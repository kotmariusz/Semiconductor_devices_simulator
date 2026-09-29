/*
 * semiconductor_3d.c  v7.3 - 3D drift-diffusion on non-uniform tensor grids
 *
 * Transport core:
 *   Poisson     : damped Newton with backtracking line search; each linear
 *                 step by PCG preconditioned with a symmetric aggregation
 *                 multigrid V-cycle (chunked D-ILU smoother, dense coarsest).
 *   Continuity  : finite-volume Scharfetter-Gummel for n and p with band-
 *                 modified potentials, solved scaled (u = c/s, row scale
 *                 1/(a_kk s_k)) so densities spanning 40+ decades get uniform
 *                 RELATIVE accuracy; right-preconditioned BiCGSTAB with a
 *                 Petrov-Galerkin aggregation multigrid (Slotboom-weighted
 *                 prolongation), ILU-only fallback.  A block-balance step
 *                 pins floating majority regions (thyristor n- base, ...).
 *   Coupling    : Gummel map on the full state (psi, phi_n, phi_p) with
 *                 Anderson acceleration (depth <= 12, RAM-aware), inexact
 *                 inner solves, quasi-Fermi potentials clamped to the span of
 *                 the terminal voltages (maximum principle).
 *   Mobility    : Arora low-field + Caughey-Thomas velocity saturation per
 *                 edge, driven by the harmonic mean of the quasi-Fermi and
 *                 electrostatic potential drops (zero at equilibrium).
 *   Continuation: adaptive steps, linear extrapolation, region-shift
 *                 predictor on the first step, MOS gates ramped first.
 *   Boundaries  : true finite-volume half cells (zero-flux Neumann) on free
 *                 surfaces; harmonic-mean face permittivity.
 *   Heterojn.   : Anderson rule (dEc = -d(chi)); per-material effective
 *                 densities enter Poisson and the SG band potentials.
 *   Gates       : oxide Robin condition C_ox (V_G' - phi_s) on a face, or a
 *                 buried metal electrode inside an insulator region.
 *   Currents    : terminal current = flux balance of the contact nodes (sums
 *                 to zero over all terminals); each terminal reports the
 *                 better resolved of its own flux and minus the sum of the
 *                 others, with a round-off resolution floor.
 *
 * Gauge: phi = 0 at the intrinsic level of the reference material (mats[0],
 * or the first semiconductor added). Applied voltages are quasi-Fermi
 * potentials: at an ohmic terminal phi_n = phi_p = V.
 * Build: gcc -O3 -march=native -ffast-math -fno-finite-math-only -fopenmp
 *        -shared -fPIC -o semiconductor_core.so semiconductor_3d.c -lm
 * (-fno-finite-math-only keeps the NaN guards of the solver alive.)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
static inline double wtime(void){ return omp_get_wtime(); }
#else
#include <time.h>
static inline double wtime(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
#endif

#define Q      1.602176634e-19
#define KB     1.380649e-23
#define EPS0   8.8541878128e-12
#define EPS_OX (3.9*EPS0)
#define HBAR   1.054571817e-34
#define M0     9.1093837015e-31
#define HPLANCK 6.62607015e-34

#define NM_MAX  12
#define NC_MAX  32
/* Carrier exponent clamp: n,p in [1e-100, 1e65] m^-3. The floor is far below
   any physically relevant density (ni of AlGaN is ~1e-8 m^-3) and keeps every
   ratio used by the scaled solvers finite. */
#define LN_MIN  (-230.0)
#define LN_MAX  150.0
#define NFLOOR  1e-100
#define PAR_MIN 16384     /* below this many nodes loops run serially */

/* size_t index arithmetic: no 32-bit overflow on very large grids */
#define IDX(i,j,l,Nx,Ny)  ((size_t)(l)*(size_t)(Nx)*(size_t)(Ny)+(size_t)(j)*(size_t)(Nx)+(size_t)(i))

enum { K_SEMI=0, K_INS=1, K_METAL=2 };

/* face-mobility storage: double by default. float halves the memory of these
   6 arrays but its 6e-8 rounding, re-applied every Gummel iteration in
   velocity-saturated channels, leaves a ~1e-8 V noise floor on phi_n. */
#define AA_MAX 12          /* max Anderson depth (memory: 24 B/node per level) */
#ifdef SEMISIM_FLOAT_MU
typedef float mu_t;
#else
typedef double mu_t;
#endif
typedef struct {
    int    id, insulator;
    double T,VT,Eg,chi,eps_r,Nc,Nv,ni;
    double mu_n_max,mu_p_max,mu_n_min,mu_p_min;
    double Nref_n,Nref_p,alpha_n,alpha_p;
    double vsat_n,vsat_p,beta_n,beta_p;
    double tau_n0,tau_p0,Cn,Cp,Brad;
    double tau_Nref;              /* SRH lifetime doping dependence (0: off)  */
    double bgn_E,bgn_N,bgn_C;     /* band-gap narrowing (Klaassen form; 0: off) */
    /* quantum confinement (MLDA): confinement masses (m0) and weights of up
       to two carrier populations (e.g. Si Delta-2 / Delta-4 valleys)      */
    double mq_e1,mq_e2,wq_e1, mq_h1,mq_h2,wq_h1;
    /* gate tunnelling: supply mass of the semiconductor side (m0); for an
       insulator the tunnelling masses of electrons / holes in it (m0)     */
    double ms_e,ms_h, mt_e,mt_h;
    /* Lombardi surface mobility (cm units, Sentaurus form); lb_on = 0: off */
    int    lb_on;
    double lbB_n,lbC_n,lbL_n,lbD_n, lbB_p,lbC_p,lbL_p,lbD_p;
    /* Derived at solve time relative to the reference material:
         n = exp(lnc + (phi - phi_n)/VT),  p = exp(lnv + (phi_p - phi)/VT)
       For the reference material lnc = lnv = ln(ni). */
    double lnc, lnv;
} Mat;

#define TUN_MAXL 6              /* max layers of a tunnelling path          */
typedef struct {
    int    face,i0,i1,j0,j1,bc;   /* face 0..5 = domain faces; 6 = volume box */
    int    k0,k1;                 /* third index range, box contacts only     */
    double V,phi_m;
    double tox, qox;              /* gate oxide thickness (m), fixed charge (C/m^2) */
    double Vprev;                 /* terminal voltage of the last solution    */
    double I, A;                  /* terminal current (A, into device), area (m^2) */
    double In, Iflr;              /* electron component; resolution floor (A) */
    double Iraw;                  /* direct contact flux (I may be the Kirchhoff complement) */
    /* gate stack of a face (Robin) gate for tunnelling: layers from the
       semiconductor surface to the metal (material index, thickness m) */
    int    snl, smat[TUN_MAXL]; double st[TUN_MAXL];
    int    ntp;                   /* tunnelling paths ending on this contact  */
} Con;

typedef struct { double *c,*xm,*xp,*ym,*yp,*zm,*zp; } Sten;

/* One multigrid level: 7-point stencil on an (Nx,Ny,Nz) tensor grid.
   Level 0 aliases the device matrix; coarse levels own their arrays. */
#define MAXLEV 24
typedef struct {
    int Nx,Ny,Nz; size_t N;
    int cx,cy,cz;                 /* coarsened along x/y/z towards next level */
    Sten A; double *ilu;
    double *x,*b,*r;
    double *pw;                   /* prolongation weight of each node to its parent */
    double *ph;                   /* reference quasi-Fermi level / VT (coarse)       */
    double *LU; int *piv; int own; /* coarsest level: dense LU */
    /* graph (CSR) levels, unstructured mode: pattern with the reverse entry of
       every entry, off-diagonal values, diagonal; towards the next level the
       aggregate of every node, its member lists and the entry map */
    int *rp,*cj,*rv; int nnz;
    double *va,*dg;
    int *agg,*mrp,*mem,*emap;
    double *wg;                   /* geometric coupling of every entry (aggregation strength) */
} Lev;

typedef struct {
    int Nx,Ny,Nz;
    double *xs,*ys,*zs;
    double *Lx,*Ly,*Lz;           /* finite-volume lengths (half cells on faces) */

    int           *mat_idx;
    unsigned char *kind;          /* K_SEMI / K_INS / K_METAL                  */
    int           *con_mask;      /* Dirichlet terminal id+1 (ohmic/Schottky/metal) */
    int           *gate_mask;     /* Robin gate id+1 (face gate on semiconductor) */
    Mat    mats[NM_MAX]; int n_mats;
    double *Nd,*Na,*Nf,*eps;      /* Nf: fixed charge density (m^-3, signed)  */
    double *bgn;                  /* band-gap narrowing, half gap / VT         */
    /* quantum confinement: ln of the MLDA density factor (<= 0), electrons /
       holes; qd = distance (m) of the node to the nearest semiconductor/
       insulator interface (1 m = none within reach)                       */
    double *qcn,*qcp,*qd;
    double *tk,*tg;               /* tunnelling sink of the carrier being solved:
                                     rate = tk*c - tg (m^3/s, 1/s); after a
                                     solve tg holds the display J_tun (A/m^2) */
    double *phi,*phi_n,*phi_p,*n,*p;
    double *Jx,*Jy,*Jz,*mu_n,*mu_p,*R_tot;
    mu_t   *fmu[6];               /* face mobility (k,k+e_ax): [0-2] n x,y,z; [3-5] p */

    /* solver work space */
    Sten   A;
    double *ilu, *wk[7];
    double *h1[3], *h2[3];        /* continuation: last converged / previous */
    double *sc, *rhs, *sol;
    double *dsc;                  /* continuity row scale D_a = diag*s_a      */
    double *xsave;                /* Krylov initial-guess backup              */
    /* Anderson history of the full Gummel state x=(psi,phi_n,phi_p), 3N long:
       differences in float (relative precision is all AA needs), the last
       residual/image in double; Gram matrix of the differences cached. */
    float  *aa_dF[AA_MAX], *aa_dG[AA_MAX];
    double *aa_f, *aa_g, aa_gram[AA_MAX][AA_MAX];
    int    aa_m, aa_n, aa_pos, aa_have;   /* allocated depth, filled, next slot */
    const double *mg_w;           /* MG level-0 restriction weights (or NULL) */
    int    mg_qsign;              /* 0: plain aggregation; +1 electrons; -1 holes */

    Lev    lev[MAXLEV]; int nlev; /* multigrid hierarchy (level 0 = A) */

    /* Robin (oxide) gate entries: node, boundary area, contact id */
    size_t *rob_k; double *rob_A; int *rob_c; int n_rob;

    Con    cons[NC_MAX]; int n_cons;
    /* gate tunnelling paths (semiconductor node -> gate metal): node, area,
       gate contact, direction towards the gate, layers (count, material,
       thickness), injection rates G_e,G_h (1/s); the supply column of each
       path (nodes inward from the interface: node, cell length, spacing to
       the next, sink coefficients K_e,K_h in m^3/s: sink = sum K c - G)    */
    size_t *tp_k; double *tp_A; int *tp_c, *tp_dir, *tp_nl, *tp_m; double *tp_t; int n_tp, cap_tp;
    double *tp_G; int *tp_c0, *tp_cn;
    size_t *col_k; double *col_l, *col_h, *col_K; int n_col, cap_col;
    double Cref;                  /* intrinsic work function of reference (eV) */
    int    ref_mat;

    int    max_iter_p,max_iter_c,max_gummel;
    double tol_p,tol_c,omega_p,omega_c;
    int    converged,conv_iter; double conv_res;
    int    allocated, have_solution, dirty, verbose;
    double vlo, vhi;              /* carrier-terminal voltage span at current stage */
    long   pcg_its, bicg_its, newton_its, gummel_its; /* diagnostics */
    double t_pois, t_cont, t_mob, t_cur, t_aa, t_mg;   /* wall-time profile (s) */
    double pcg_rtol, ntol_fac, bicg_rtol;   /* inexact inner-solve controls */
    int    nthreads;
    int    models;                /* bit0 BGN, bit1 tau(N), bit2 v-sat, bit3 MLDA
                                     quantum confinement, bit4 Lombardi surface
                                     mobility, bit5 gate tunnelling          */
    /* unstructured (graph) mesh, unstr = 1: nodes 0..N-1 (Nx = N, Ny = Nz =
       1); neighbours from a symmetric CSR graph carrying the box-method
       geometry: per directed edge the neighbour, reverse edge, undirected id,
       length (m) and dual face area (m^2); per node the control volume (m^3)
       and position (m).  Matrix off-diagonals live per directed edge (gval),
       edge mobilities per undirected edge (gmu). */
    int    unstr;
    int    *grp, *gcol, *grev, *gue; int nE, nUE;
    double *gh, *ga, *gvol, *gpos, *gval;
    mu_t   *gmu[2];
    double *glsq;                 /* 6 per node: inverse least-squares gradient matrix */
    double *ggrad;                /* 3 per node: potential gradient (work)             */
    double *gwall;                /* 6 per node: two interface walls (distance, cell range a..b) */
    int    cnn[NC_MAX]; int *cnodes[NC_MAX]; double *careas[NC_MAX];  /* contact node lists */
} Dev;

static Dev gdev;

/* ------------------------------------------------------------ materials */
static void mat_init(Mat *m, int id, double T){
    memset(m,0,sizeof(*m));
    m->id=id; m->T=T; m->VT=KB*T/Q; m->insulator=0;
    switch(id){
    case 0: /* Si */
        m->Eg=1.170-4.73e-4*T*T/(T+636.); m->chi=4.05; m->eps_r=11.7;
        m->Nc=2.86e25*pow(T/300.,1.58); m->Nv=3.10e25*pow(T/300.,1.85);
        m->mu_n_max=0.1417; m->mu_p_max=0.0470; m->mu_n_min=0.0088; m->mu_p_min=0.00543;
        m->Nref_n=1.30e23; m->Nref_p=2.35e23; m->alpha_n=0.88; m->alpha_p=0.88;
        m->vsat_n=1.07e5; m->vsat_p=8.37e4; m->beta_n=1.109; m->beta_p=1.213;
        m->tau_n0=1e-5; m->tau_p0=1e-5; m->Cn=2.8e-43; m->Cp=9.9e-44; m->Brad=1.8e-21;
        m->tau_Nref=1e22;                           /* Scharfetter, 1e16 cm^-3 */
        m->bgn_E=6.92e-3; m->bgn_N=1.3e23; m->bgn_C=0.5;  /* Klaassen 1992 */
        /* MLDA, (100) interface: every valley keeps its Boltzmann share of
           the density of states and its own confinement mass - electrons in
           the 2 Delta-2 valleys (m_l = 0.916, 1/3) and the 4 Delta-4 valleys
           (m_t = 0.19, 2/3); holes heavy (m_z ~0.29, 0.84) and light (~0.20) */
        m->mq_e1=0.916; m->mq_e2=0.19; m->wq_e1=1./3.; m->mq_h1=0.29; m->mq_h2=0.20; m->wq_h1=0.84;
        m->ms_e=0.50; m->ms_h=0.50;
        /* Lombardi (1988) as parameterised in Sentaurus: B/F + C N^l / F^(1/3); delta/F^2 */
        m->lb_on=1;
        m->lbB_n=4.75e7; m->lbC_n=580.;  m->lbL_n=0.125;  m->lbD_n=5.82e14;
        m->lbB_p=9.925e6; m->lbC_p=2947.; m->lbL_p=0.0317; m->lbD_p=2.0546e14;
        break;
    case 1: /* Ge */
        m->Eg=0.742-4.77e-4*T*T/(T+235.); m->chi=4.00; m->eps_r=16.2;
        m->Nc=1.04e25*pow(T/300.,1.5); m->Nv=6.0e25*pow(T/300.,1.5);
        m->mu_n_max=0.39; m->mu_p_max=0.19; m->mu_n_min=0.040; m->mu_p_min=0.020;
        m->Nref_n=1.30e23; m->Nref_p=1.00e23; m->alpha_n=0.56; m->alpha_p=0.43;
        m->vsat_n=6e4; m->vsat_p=5.4e4; m->beta_n=1.; m->beta_p=1.;
        m->tau_n0=1e-3; m->tau_p0=1e-3; m->Cn=1e-42; m->Cp=1e-42; m->Brad=6.4e-20; break;
    case 2: /* GaAs */
        m->Eg=1.519-5.405e-4*T*T/(T+204.); m->chi=4.07; m->eps_r=12.9;
        m->Nc=4.37e23*pow(T/300.,1.5); m->Nv=8.68e24*pow(T/300.,1.5);
        m->mu_n_max=0.80; m->mu_p_max=0.040; m->mu_n_min=0.050; m->mu_p_min=0.0020;
        m->Nref_n=1.70e23; m->Nref_p=2.75e23; m->alpha_n=0.394; m->alpha_p=0.394;
        m->vsat_n=7.7e4; m->vsat_p=7.7e4; m->beta_n=1.; m->beta_p=1.;
        m->tau_n0=1e-9; m->tau_p0=1e-9; m->Cn=1e-43; m->Cp=1e-43; m->Brad=1.5e-16; break;
    case 3: /* InP */
        m->Eg=1.344-4.5e-4*T*T/(T+327.); m->chi=4.38; m->eps_r=12.5;
        m->Nc=5.7e23*pow(T/300.,1.5); m->Nv=1.1e25*pow(T/300.,1.5);
        m->mu_n_max=0.45; m->mu_p_max=0.015; m->mu_n_min=0.040; m->mu_p_min=0.0010;
        m->Nref_n=3.00e23; m->Nref_p=4.90e23; m->alpha_n=0.40; m->alpha_p=0.40;
        m->vsat_n=1e5; m->vsat_p=5e4; m->beta_n=1.; m->beta_p=1.;
        m->tau_n0=1e-8; m->tau_p0=1e-8; m->Cn=1e-43; m->Cp=1e-43; m->Brad=8.5e-17; break;
    case 4: /* 4H-SiC */
        m->Eg=3.265-6.5e-4*T*T/(T+1300.); m->chi=3.8; m->eps_r=9.7;
        m->Nc=1.73e25*pow(T/300.,1.5); m->Nv=2.64e25*pow(T/300.,1.5);
        m->mu_n_max=0.095; m->mu_p_max=0.0125; m->mu_n_min=0.0040; m->mu_p_min=0.0015;
        m->Nref_n=1.94e23; m->Nref_p=2.35e24; m->alpha_n=0.61; m->alpha_p=0.34;
        m->vsat_n=2.2e5; m->vsat_p=1.6e5; m->beta_n=1.2; m->beta_p=1.;
        m->tau_n0=1e-6; m->tau_p0=1e-7; m->Cn=5e-44; m->Cp=2e-44; m->Brad=1.5e-23; break;
    case 5: /* GaN */
        m->Eg=3.507-9.09e-4*T*T/(T+830.); m->chi=4.1; m->eps_r=8.9;
        m->Nc=2.3e24*pow(T/300.,1.5); m->Nv=4.6e25*pow(T/300.,1.5);
        m->mu_n_max=0.120; m->mu_p_max=0.0030; m->mu_n_min=0.0055; m->mu_p_min=0.00030;
        m->Nref_n=1.00e23; m->Nref_p=3.00e23; m->alpha_n=1.0; m->alpha_p=2.0;
        m->vsat_n=2.5e5; m->vsat_p=2e5; m->beta_n=1.; m->beta_p=1.;
        m->tau_n0=1e-8; m->tau_p0=1e-8; m->Cn=1e-43; m->Cp=1e-43; m->Brad=1.1e-16; break;
    case 6: /* Al(x)Ga(1-x)As, x=0.3: Eg(300K)=1.424+1.247x=1.798 eV (direct);
               band offsets dEc:dEv = 65:35 vs GaAs -> dEc=0.243 eV, chi=3.827 */
        m->Eg=1.882-5.58e-4*T*T/(T+295.); m->chi=3.827; m->eps_r=12.24;
        m->Nc=7.2e23*pow(T/300.,1.5); m->Nv=1.05e25*pow(T/300.,1.5);
        m->mu_n_max=0.35; m->mu_p_max=0.025; m->mu_n_min=0.020; m->mu_p_min=0.0015;
        m->Nref_n=6.00e22; m->Nref_p=1.00e23; m->alpha_n=0.394; m->alpha_p=0.394;
        m->vsat_n=6e4; m->vsat_p=5e4; m->beta_n=1.; m->beta_p=1.;
        m->tau_n0=1e-9; m->tau_p0=1e-9; m->Cn=1e-43; m->Cp=1e-43; m->Brad=5e-17; break;
    case 7: /* Al(x)Ga(1-x)N, x=0.25 - HEMT barrier.
               Eg: Vegard + bowing b=1.0 eV (Vurgaftman & Meyer 2003) -> 3.92 eV @300K
               dEc = 0.7*dEg vs GaN -> chi = 4.10-0.34 = 3.76 eV */
        m->Eg=3.990-9.09e-4*T*T/(T+830.); m->chi=3.76; m->eps_r=8.8;
        m->Nc=3.2e24*pow(T/300.,1.5); m->Nv=4.6e25*pow(T/300.,1.5);
        m->mu_n_max=0.030; m->mu_p_max=0.0010; m->mu_n_min=0.0030; m->mu_p_min=0.00020;
        m->Nref_n=1.00e23; m->Nref_p=3.00e23; m->alpha_n=1.0; m->alpha_p=2.0;
        m->vsat_n=1.5e5; m->vsat_p=1e5; m->beta_n=1.; m->beta_p=1.;
        m->tau_n0=1e-9; m->tau_p0=1e-9; m->Cn=1e-43; m->Cp=1e-43; m->Brad=1e-16; break;
    /* Insulators: chi sets the band offsets to the semiconductors (SiO2/Si:
       dEc 3.10, dEv 4.78 eV); mt_* are the tunnelling masses in the layer. */
    case 8: /* SiO2 (insulator) */
        m->insulator=1; m->Eg=9.0; m->chi=0.95; m->eps_r=3.9; m->mt_e=0.42; m->mt_h=0.33; break;
    case 9: /* Al2O3 (insulator, ALD): dEc(Si) 2.7 eV */
        m->insulator=1; m->Eg=6.8; m->chi=1.35; m->eps_r=9.0; m->mt_e=0.35; m->mt_h=0.35; break;
    case 10: /* HfO2 (high-k): k 22, Eg 5.8 eV, dEc(Si) 1.5 eV, dEv(Si) ~3.2 eV */
        m->insulator=1; m->Eg=5.8; m->chi=2.55; m->eps_r=22.0; m->mt_e=0.18; m->mt_h=0.20; break;
    case 11: /* Si3N4: k 7.5, Eg 5.1 eV, dEc(Si) 2.0 eV */
        m->insulator=1; m->Eg=5.1; m->chi=2.05; m->eps_r=7.5; m->mt_e=0.50; m->mt_h=0.50; break;
    default:
        fprintf(stderr,"semisim: WARNING unknown material id %d, using Si\n",id);
        mat_init(m,0,T); return;
    }
    if(!m->insulator && !(m->mq_e1>0.)){    /* generic: DOS-like single masses */
        static const double me[8]={0.916,0.082,0.067,0.08,0.37,0.20,0.09,0.23};
        static const double mh[8]={0.29, 0.28, 0.50, 0.60,1.00,1.00,0.55,1.00};
        m->mq_e1=m->mq_e2=me[id]; m->wq_e1=1.; m->mq_h1=m->mq_h2=mh[id]; m->wq_h1=1.;
        m->ms_e=me[id]<0.5?me[id]:0.5; m->ms_h=0.5;
    }
    m->ni = m->insulator ? 0. : sqrt(m->Nc*m->Nv)*exp(-m->Eg/(2.*m->VT));
}

/* Band offsets: pick the reference (first semiconductor), then give every
   semiconductor its effective ln(Nc), ln(Nv) in the common gauge. */
static void mat_prepare(Dev *d){
    d->ref_mat=0;
    for(int i=0;i<d->n_mats;i++) if(!d->mats[i].insulator){ d->ref_mat=i; break; }
    const Mat *r=&d->mats[d->ref_mat];
    d->Cref = r->insulator ? 4.6 : r->chi + 0.5*r->Eg + 0.5*r->VT*log(r->Nc/r->Nv);
    for(int i=0;i<d->n_mats;i++){
        Mat *m=&d->mats[i];
        if(m->insulator){ m->lnc=m->lnv=-1e300; continue; }
        m->lnc = log(m->Nc) + (m->chi - d->Cref)/m->VT;
        m->lnv = log(m->Nv) + (d->Cref - m->chi - m->Eg)/m->VT;
    }
}

/* Arora / Caughey-Thomas low-field mobility + Caughey-Thomas saturation.
   E is the DRIVING field for that carrier (|grad phi_n| or |grad phi_p|). */
static double mob_low(const Mat *m, double Ntot, int hole){
    if(Ntot<1.)Ntot=1.;
    double tr=m->T/300.;
    double mumin = (hole?m->mu_p_min:m->mu_n_min)*pow(tr,-0.57);
    double mumax = (hole?m->mu_p_max:m->mu_n_max)*pow(tr,hole?-2.23:-2.33);
    double Nref  = (hole?m->Nref_p:m->Nref_n)*pow(tr,2.546);
    double alpha = (hole?m->alpha_p:m->alpha_n)*pow(tr,-0.146);
    if(mumax<mumin) mumax=mumin;
    return mumin + (mumax-mumin)/(1.+pow(Ntot/Nref,alpha));
}
static double mob_hf(double mu0, double E, double vsat, double beta){
    double r=mu0*fabs(E)/vsat;
    if(r<1e-12) return mu0;
    return mu0/pow(1.+pow(r,beta),1./beta);
}

static inline double B(double x){                 /* Bernoulli x/(e^x-1) */
    if(fabs(x)<1e-6) return 1.-x/2.+x*x/12.;
    if(x> 700.) return 0.;
    if(x<-700.) return -x;
    return x/expm1(x);
}

static inline double clampe(double e){ return e<LN_MIN?LN_MIN:(e>LN_MAX?LN_MAX:e); }
static inline double n_of(const Mat *m, double phi, double phin){
    return exp(clampe(m->lnc+(phi-phin)/m->VT));
}
static inline double p_of(const Mat *m, double phi, double phip){
    return exp(clampe(m->lnv+(phip-phi)/m->VT));
}

/* SRH (midgap trap) + Auger + radiative; derivatives for linearisation */
static inline void recomb(const Mat *m, double ni, double tf, double n, double p,
                          double *R, double *dRn, double *dRp){
    double ni2=ni*ni, np=n*p-ni2;
    double den=tf*(m->tau_p0*(n+ni)+m->tau_n0*(p+ni))+1e-300;
    double Rs=np/den;
    double dsn=(p*den-np*tf*m->tau_p0)/(den*den);
    double dsp=(n*den-np*tf*m->tau_n0)/(den*den);
    double Au=m->Cn*n+m->Cp*p;
    *R  = Rs + Au*np + m->Brad*np;
    *dRn= dsn + m->Cn*np + Au*p + m->Brad*p;
    *dRp= dsp + m->Cp*np + Au*n + m->Brad*n;
}

/* ----------------------------------------------------------- geometry */
static inline double fv_len(const double *c, int i, int N){
    if(N<2) return 1.;
    if(i==0)   return 0.5*(c[1]-c[0]);
    if(i==N-1) return 0.5*(c[N-1]-c[N-2]);
    return 0.5*(c[i+1]-c[i-1]);
}

/* Neighbour enumeration: 6 directions, returns 0 if outside the domain.
   dir: 0 x-,1 x+,2 y-,3 y+,4 z-,5 z+.  h = node spacing, area = face area. */
static inline int nbr(const Dev *d, int i, int j, int l, int dir,
                      size_t *kn, double *h, double *area){
    int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz;
    switch(dir){
    case 0: if(i==0)    return 0; *kn=IDX(i-1,j,l,Nx,Ny); *h=d->xs[i]-d->xs[i-1]; *area=d->Ly[j]*d->Lz[l]; return 1;
    case 1: if(i==Nx-1) return 0; *kn=IDX(i+1,j,l,Nx,Ny); *h=d->xs[i+1]-d->xs[i]; *area=d->Ly[j]*d->Lz[l]; return 1;
    case 2: if(j==0)    return 0; *kn=IDX(i,j-1,l,Nx,Ny); *h=d->ys[j]-d->ys[j-1]; *area=d->Lx[i]*d->Lz[l]; return 1;
    case 3: if(j==Ny-1) return 0; *kn=IDX(i,j+1,l,Nx,Ny); *h=d->ys[j+1]-d->ys[j]; *area=d->Lx[i]*d->Lz[l]; return 1;
    case 4: if(l==0)    return 0; *kn=IDX(i,j,l-1,Nx,Ny); *h=d->zs[l]-d->zs[l-1]; *area=d->Lx[i]*d->Ly[j]; return 1;
    default:if(l==Nz-1) return 0; *kn=IDX(i,j,l+1,Nx,Ny); *h=d->zs[l+1]-d->zs[l]; *area=d->Lx[i]*d->Ly[j]; return 1;
    }
}
static inline double volume(const Dev *d, int i, int j, int l){
    return d->Lx[i]*d->Ly[j]*d->Lz[l];
}
static inline void unidx(const Dev *d, size_t k, int *i, int *j, int *l){
    size_t NxNy=(size_t)d->Nx*(size_t)d->Ny;
    *l=(int)(k/NxNy); size_t r=k-(size_t)(*l)*NxNy;
    *j=(int)(r/(size_t)d->Nx); *i=(int)(r-(size_t)(*j)*(size_t)d->Nx);
}

/* ------------------------------------------------ node-level band data
   Effective densities include band-gap narrowing, split symmetrically
   (Ec down, Ev up by dEg/2): ln Nc_eff = lnc + bgn, ln Nv_eff = lnv + bgn. */
static inline const Mat *MK(const Dev *d, size_t k){ return &d->mats[d->mat_idx[k]]; }
/* ln of the effective densities of states at node k: material + band-gap
   narrowing + quantum-confinement factor (electrons qcn, holes qcp) */
static inline double lncK(const Dev *d, size_t k){ return MK(d,k)->lnc+d->bgn[k]+d->qcn[k]; }
static inline double lnvK(const Dev *d, size_t k){ return MK(d,k)->lnv+d->bgn[k]+d->qcp[k]; }
static inline double niK (const Dev *d, size_t k){ return MK(d,k)->ni*exp(d->bgn[k]+0.5*(d->qcn[k]+d->qcp[k])); }
static inline double nK(const Dev *d, size_t k, double psi, double qn){
    const Mat *m=MK(d,k); return exp(clampe(m->lnc+d->bgn[k]+d->qcn[k]+(psi-qn)/m->VT));
}
static inline double pK(const Dev *d, size_t k, double psi, double qp){
    const Mat *m=MK(d,k); return exp(clampe(m->lnv+d->bgn[k]+d->qcp[k]+(qp-psi)/m->VT));
}
/* band-potential differences of an edge a->b (units of VT): electrons
   D = dpsi/VT + d(ln Nc_eff), holes D = dpsi/VT - d(ln Nv_eff) */
static inline double dlnc(const Dev *d, size_t a, size_t b){ return lncK(d,b)-lncK(d,a); }
static inline double dlnv(const Dev *d, size_t a, size_t b){ return lnvK(d,b)-lnvK(d,a); }
static inline double tauf(const Dev *d, size_t k){
    const Mat *m=MK(d,k);
    if(!(d->models&2) || !(m->tau_Nref>0.)) return 1.;
    return 1./(1.+(d->Nd[k]+d->Na[k])/m->tau_Nref);
}
static inline void recombK(const Dev *d, size_t k, double *R, double *dRn, double *dRp){
    recomb(MK(d,k),niK(d,k),tauf(d,k),d->n[k],d->p[k],R,dRn,dRp);
}
static double bgn_half(const Mat *m, double N){       /* dEg/(2 VT) */
    if(!(m->bgn_E>0.) || N<1.) return 0.;
    double L=log(N/m->bgn_N);
    return 0.5*m->bgn_E*(L+sqrt(L*L+m->bgn_C))/m->VT;
}

/* ----------------------------------------------------------- allocation */
#define NODE_D(d) { &d->Nd,&d->Na,&d->Nf,&d->eps,&d->bgn,&d->qcn,&d->qcp,&d->qd,&d->tk,&d->tg, \
    &d->phi,&d->phi_n,&d->phi_p,&d->n,&d->p, \
    &d->Jx,&d->Jy,&d->Jz,&d->mu_n,&d->mu_p,&d->R_tot, \
    &d->A.c,&d->A.xm,&d->A.xp,&d->A.ym,&d->A.yp,&d->A.zm,&d->A.zp,&d->ilu, \
    &d->wk[0],&d->wk[1],&d->wk[2],&d->wk[3],&d->wk[4],&d->wk[5],&d->wk[6], \
    &d->h1[0],&d->h1[1],&d->h1[2],&d->h2[0],&d->h2[1],&d->h2[2], \
    &d->sc,&d->rhs,&d->sol,&d->dsc,&d->xsave }
#define N_NODE_D 47

static void mg_free(Dev *d);
static int  mg_alloc(Dev *d);
static void tun_free(Dev *d){
    free(d->tp_k); free(d->tp_A); free(d->tp_c); free(d->tp_dir); free(d->tp_nl); free(d->tp_m);
    free(d->tp_t); free(d->tp_G); free(d->tp_c0); free(d->tp_cn);
    free(d->col_k); free(d->col_l); free(d->col_h); free(d->col_K);
    d->tp_k=NULL; d->tp_A=d->tp_t=d->tp_G=NULL; d->tp_c=d->tp_dir=d->tp_nl=d->tp_m=d->tp_c0=d->tp_cn=NULL;
    d->col_k=NULL; d->col_l=d->col_h=d->col_K=NULL;
    d->n_tp=d->cap_tp=d->n_col=d->cap_col=0;
}

static void aa_free(Dev *d){
    for(int i=0;i<AA_MAX;i++){ free(d->aa_dF[i]); free(d->aa_dG[i]); d->aa_dF[i]=d->aa_dG[i]=NULL; }
    free(d->aa_f); free(d->aa_g); d->aa_f=d->aa_g=NULL; d->aa_m=0;
}
/* Anderson history for a 3N state; degrade the depth if memory is short. */
static int aa_alloc(Dev *d, size_t N, int want){
    aa_free(d);
    size_t L=3*N;
    d->aa_f=malloc(L*sizeof(double)); d->aa_g=malloc(L*sizeof(double));
    if(!d->aa_f||!d->aa_g){ aa_free(d); return 0; }
    for(int m=want; m>=2; m=(m>4? m-4 : m-2)){
        int ok=1;
        for(int i=0;i<m;i++){
            d->aa_dF[i]=malloc(L*sizeof(float)); d->aa_dG[i]=malloc(L*sizeof(float));
            if(!d->aa_dF[i]||!d->aa_dG[i]){ ok=0; break; }
        }
        if(ok){ d->aa_m=m; return 1; }
        for(int i=0;i<AA_MAX;i++){ free(d->aa_dF[i]); free(d->aa_dG[i]); d->aa_dF[i]=d->aa_dG[i]=NULL; }
    }
    d->aa_m=0; return 1;            /* plain Gummel still works */
}

static void gfree(Dev *d){
    free(d->grp); free(d->gcol); free(d->grev); free(d->gue); free(d->gh); free(d->ga);
    free(d->gvol); free(d->gpos); free(d->gval); free(d->gmu[0]); free(d->gmu[1]);
    free(d->glsq); free(d->ggrad); free(d->gwall);
    d->grp=d->gcol=d->grev=d->gue=NULL; d->gh=d->ga=d->gvol=d->gpos=d->gval=NULL;
    d->gmu[0]=d->gmu[1]=NULL; d->glsq=d->ggrad=d->gwall=NULL;
    for(int c=0;c<NC_MAX;c++){ free(d->cnodes[c]); free(d->careas[c]); d->cnodes[c]=NULL; d->careas[c]=NULL; d->cnn[c]=0; }
    d->unstr=0; d->nE=d->nUE=0;
}

static void free_dev(Dev *d){
    mg_free(d);
    aa_free(d);
    gfree(d);
    double **pp[N_NODE_D]=NODE_D(d);
    for(int q=0;q<N_NODE_D;q++){ free(*pp[q]); *pp[q]=NULL; }
    free(d->xs); free(d->ys); free(d->zs); free(d->Lx); free(d->Ly); free(d->Lz);
    free(d->mat_idx); free(d->kind); free(d->con_mask); free(d->gate_mask);
    free(d->rob_k); free(d->rob_A); free(d->rob_c);
    tun_free(d);
    for(int q=0;q<6;q++){ free(d->fmu[q]); d->fmu[q]=NULL; }
    d->xs=d->ys=d->zs=d->Lx=d->Ly=d->Lz=NULL;
    d->mat_idx=d->con_mask=d->gate_mask=NULL; d->kind=NULL;
    d->rob_k=NULL; d->rob_A=NULL; d->rob_c=NULL; d->n_rob=0;
    d->allocated=0; d->have_solution=0;
}

static int alloc_dev(Dev *d, int Nx, int Ny, int Nz){
    free_dev(d);
    if(Nx<2||Ny<2||Nz<2) return 0;
    size_t N=(size_t)Nx*(size_t)Ny*(size_t)Nz;
    if(N>2000000000ULL) return 0;
    int ok=1;
    double **pp[N_NODE_D]=NODE_D(d);
    for(int q=0;q<N_NODE_D;q++){ *pp[q]=calloc(N,sizeof(double)); if(!*pp[q]) ok=0; }
    d->xs=calloc(Nx,sizeof(double)); d->Lx=calloc(Nx,sizeof(double));
    d->ys=calloc(Ny,sizeof(double)); d->Ly=calloc(Ny,sizeof(double));
    d->zs=calloc(Nz,sizeof(double)); d->Lz=calloc(Nz,sizeof(double));
    d->mat_idx=calloc(N,sizeof(int)); d->con_mask=calloc(N,sizeof(int));
    d->gate_mask=calloc(N,sizeof(int)); d->kind=calloc(N,1);
    for(int q=0;q<6;q++){ d->fmu[q]=calloc(N,sizeof(mu_t)); if(!d->fmu[q]) ok=0; }
    if(!d->xs||!d->ys||!d->zs||!d->Lx||!d->Ly||!d->Lz||!d->mat_idx||!d->con_mask
       ||!d->gate_mask||!d->kind) ok=0;
    if(!ok){ free_dev(d); return 0; }
    d->Nx=Nx; d->Ny=Ny; d->Nz=Nz; d->allocated=1; d->have_solution=0; d->dirty=1;
    for(size_t k=0;k<N;k++) d->qd[k]=1.;     /* no interface nearby */
    if(!mg_alloc(d)){ free_dev(d); return 0; }
    /* Anderson depth: AA_MAX unless physical RAM is short (the history is an
       accelerator only: 24 B/node per level) */
    int want=AA_MAX;
#ifdef _SC_PHYS_PAGES
    { long pg=sysconf(_SC_PHYS_PAGES), ps=sysconf(_SC_PAGE_SIZE);
      if(pg>0 && ps>0){
          double ram=(double)pg*(double)ps;
          double base=(double)N*(N_NODE_D*8.+6.*sizeof(mu_t)+13.+16.+48.);
          int fit=(int)floor((0.75*ram-base)/((double)N*24.));
          if(fit<want) want=fit;
      } }
#endif
    if(!aa_alloc(d,N,want)){ free_dev(d); return 0; }
    return 1;
}

static inline size_t NN(const Dev *d){ return (size_t)d->Nx*(size_t)d->Ny*(size_t)d->Nz; }

static int gmg_alloc(Dev *d);
static void g_lsq_setup(Dev *d);
/* Unstructured device: N nodes, E directed edges (symmetric CSR: rp, cj),
   per directed edge the length h (m) and the dual face area ar (m^2), per
   node the control volume (m^3) and position pos (3N, m).  Nx = N, Ny = Nz
   = 1, so every node-wise loop of the tensor code still applies. */
static int alloc_dev_graph(Dev *d, int N, int E, const int *rp, const int *cj, const double *h,
                           const double *ar, const double *vol, const double *pos){
    free_dev(d);
    if(N<2||E<1||rp[0]!=0||rp[N]!=E) return 0;
    int ok=1;
    double **pp[N_NODE_D]=NODE_D(d);
    for(int q=0;q<N_NODE_D;q++){ *pp[q]=calloc((size_t)N,sizeof(double)); if(!*pp[q]) ok=0; }
    d->xs=calloc((size_t)N,sizeof(double)); d->Lx=calloc((size_t)N,sizeof(double));
    d->ys=calloc(1,sizeof(double)); d->Ly=calloc(1,sizeof(double));
    d->zs=calloc(1,sizeof(double)); d->Lz=calloc(1,sizeof(double));
    d->mat_idx=calloc((size_t)N,sizeof(int)); d->con_mask=calloc((size_t)N,sizeof(int));
    d->gate_mask=calloc((size_t)N,sizeof(int)); d->kind=calloc((size_t)N,1);
    for(int q=0;q<6;q++){ d->fmu[q]=calloc(1,sizeof(mu_t)); if(!d->fmu[q]) ok=0; }
    d->grp=malloc((size_t)(N+1)*sizeof(int)); d->gcol=malloc((size_t)E*sizeof(int));
    d->grev=malloc((size_t)E*sizeof(int)); d->gue=malloc((size_t)E*sizeof(int));
    d->gh=malloc((size_t)E*8); d->ga=malloc((size_t)E*8); d->gval=calloc((size_t)E,8);
    d->gvol=malloc((size_t)N*8); d->gpos=malloc((size_t)3*N*8);
    d->glsq=calloc((size_t)6*N,8); d->ggrad=calloc((size_t)3*N,8); d->gwall=malloc((size_t)6*N*8);
    if(!d->xs||!d->ys||!d->zs||!d->Lx||!d->Ly||!d->Lz||!d->mat_idx||!d->con_mask||!d->gate_mask||!d->kind
       ||!d->grp||!d->gcol||!d->grev||!d->gue||!d->gh||!d->ga||!d->gval||!d->gvol||!d->gpos
       ||!d->glsq||!d->ggrad||!d->gwall) ok=0;
    if(!ok){ free_dev(d); return 0; }
    memcpy(d->grp,rp,(size_t)(N+1)*sizeof(int)); memcpy(d->gcol,cj,(size_t)E*sizeof(int));
    memcpy(d->gh,h,(size_t)E*8); memcpy(d->ga,ar,(size_t)E*8);
    memcpy(d->gvol,vol,(size_t)N*8); memcpy(d->gpos,pos,(size_t)3*N*8);
    for(int k=0;k<N;k++){
        if(rp[k+1]<rp[k] || !(vol[k]>0.)) ok=0;
        for(int e=rp[k];e<rp[k+1] && ok;e++){
            if(cj[e]<0||cj[e]>=N||cj[e]==k||!(h[e]>0.)||!(ar[e]>=0.)) ok=0;
        }
    }
    /* reverse edges and undirected ids */
    int nue=0;
    for(int k=0;k<N && ok;k++) for(int e=rp[k];e<rp[k+1];e++){
        int j=cj[e], r=-1;
        for(int t=rp[j];t<rp[j+1];t++) if(cj[t]==k){ r=t; break; }
        if(r<0){ ok=0; break; }
        d->grev[e]=r;
        if(k<j) d->gue[e]=nue++;
    }
    if(ok) for(int k=0;k<N;k++) for(int e=rp[k];e<rp[k+1];e++) if(cj[e]<k) d->gue[e]=d->gue[d->grev[e]];
    if(!ok){ free_dev(d); return 0; }
    d->nE=E; d->nUE=nue;
    d->gmu[0]=calloc((size_t)(nue>0?nue:1),sizeof(mu_t)); d->gmu[1]=calloc((size_t)(nue>0?nue:1),sizeof(mu_t));
    if(!d->gmu[0]||!d->gmu[1]){ free_dev(d); return 0; }
    for(int k=0;k<N;k++){ d->xs[k]=pos[3*k]; d->Lx[k]=vol[k]; }
    d->ys[0]=d->zs[0]=0.; d->Ly[0]=d->Lz[0]=1.;
    d->Nx=N; d->Ny=1; d->Nz=1; d->unstr=1; d->allocated=1; d->have_solution=0; d->dirty=1;
    for(int k=0;k<N;k++){ d->qd[k]=1.; for(int q=0;q<6;q++) d->gwall[6*k+q]=1.; }
    if(!gmg_alloc(d)){ free_dev(d); return 0; }
    int want=AA_MAX;
#ifdef _SC_PHYS_PAGES
    { long pg=sysconf(_SC_PHYS_PAGES), ps=sysconf(_SC_PAGE_SIZE);
      if(pg>0 && ps>0){
          double ram=(double)pg*(double)ps;
          double base=(double)N*(N_NODE_D*8.+13.+16.+48.+120.)+(double)E*40.;
          int fit=(int)floor((0.75*ram-base)/((double)N*24.));
          if(fit<want) want=fit;
      } }
#endif
    if(!aa_alloc(d,(size_t)N,want)){ free_dev(d); return 0; }
    return 1;
}

/* ------------------------------------------------------- device setup */
/* Iterate over the nodes of contact c: calls body with (i,j,l,k, boundary area).
   Face contacts use the v6 index convention; face 6 is a volume box. */
#define FOR_CON_NODES(d,cn,...) do{ \
    int _Nx=(d)->Nx,_Ny=(d)->Ny,_Nz=(d)->Nz; \
    int _ia=0,_ib=_Nx-1,_ja=0,_jb=_Ny-1,_la=0,_lb=_Nz-1; \
    switch((cn)->face){ \
    case 0: _ia=_ib=0;       _ja=(cn)->i0;_jb=(cn)->i1;_la=(cn)->j0;_lb=(cn)->j1; break; \
    case 1: _ia=_ib=_Nx-1;   _ja=(cn)->i0;_jb=(cn)->i1;_la=(cn)->j0;_lb=(cn)->j1; break; \
    case 2: _ja=_jb=0;       _ia=(cn)->i0;_ib=(cn)->i1;_la=(cn)->j0;_lb=(cn)->j1; break; \
    case 3: _ja=_jb=_Ny-1;   _ia=(cn)->i0;_ib=(cn)->i1;_la=(cn)->j0;_lb=(cn)->j1; break; \
    case 4: _la=_lb=0;       _ia=(cn)->i0;_ib=(cn)->i1;_ja=(cn)->j0;_jb=(cn)->j1; break; \
    case 5: _la=_lb=_Nz-1;   _ia=(cn)->i0;_ib=(cn)->i1;_ja=(cn)->j0;_jb=(cn)->j1; break; \
    default: _ia=(cn)->i0;_ib=(cn)->i1;_ja=(cn)->j0;_jb=(cn)->j1;_la=(cn)->k0;_lb=(cn)->k1; \
    } \
    if(_ia<0)_ia=0; if(_ja<0)_ja=0; if(_la<0)_la=0; \
    if(_ib>_Nx-1)_ib=_Nx-1; if(_jb>_Ny-1)_jb=_Ny-1; if(_lb>_Nz-1)_lb=_Nz-1; \
    for(int l=_la;l<=_lb;l++) for(int j=_ja;j<=_jb;j++) for(int i=_ia;i<=_ib;i++){ \
        size_t k=IDX(i,j,l,_Nx,_Ny); \
        double area = ((cn)->face<=1) ? (d)->Ly[j]*(d)->Lz[l] : \
                      ((cn)->face<=3) ? (d)->Lx[i]*(d)->Lz[l] : \
                      ((cn)->face<=5) ? (d)->Lx[i]*(d)->Ly[j] : 0.; \
        (void)area; __VA_ARGS__ \
    } }while(0)

/* MLDA (modified local-density approximation, Paasch & Uebensee 1982) for
   Boltzmann carriers: in front of a hard wall the density of a population
   with confinement mass m falls off as 1 - exp(-(z/lambda)^2),
   lambda = hbar/sqrt(2 m kT) (1.27 nm for Si electrons, m_l = 0.916).
   The factor is averaged over the node's cell, z in [a,b] from the wall;
   up to two populations (weights w, 1-w). */
static double mlda_pop(double a, double b, double lam){
    double u=a/lam, v=b/lam;
    if(v<0.05){                                   /* series: no cancellation */
        double d=v-u; if(d<1e-12) return u*u;
        return ((v*v*v-u*u*u)/3.-(pow(v,5)-pow(u,5))/10.+(pow(v,7)-pow(u,7))/42.)/d;
    }
    if(v-u<1e-9*v) return 1.-exp(-u*u);
    return 1.-0.5*sqrt(M_PI)*(erf(v)-erf(u))/(v-u);
}
static double mlda_factor(const Mat *m, double a, double b, int hole){
    double kT=KB*m->T, m1=hole?m->mq_h1:m->mq_e1, m2=hole?m->mq_h2:m->mq_e2, w=hole?m->wq_h1:m->wq_e1;
    if(!(m1>0.)) return 1.;
    if(a<0.) a=0.; if(b<a) b=a;
    double f=w*mlda_pop(a,b,HBAR/sqrt(2.*m1*M0*kT));
    if(w<1. && m2>0.) f+=(1.-w)*mlda_pop(a,b,HBAR/sqrt(2.*m2*M0*kT));
    return f;
}

/* Interfaces on the tensor mesh.  Along each axis every semiconductor node
   finds the nearest wall on either side within its contiguous run of
   semiconductor nodes: a semiconductor-insulator node pair (wall at the
   midpoint - the box-method interface) or a Robin gate on the domain face
   (wall at the node).  qd = distance to the nearest wall (Lombardi); with
   MLDA on, qcn/qcp = sum over the walls of ln(cell-averaged factor), i.e.
   the product of the single-wall factors (exact for a rectangular corner
   with Boltzmann statistics). */
#define WALL_REACH 40e-9                          /* factor = 1 beyond this */
static void walls_tensor(Dev *d, int qc){
    const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz; size_t N=NN(d);
    for(size_t k=0;k<N;k++){ d->qd[k]=1.; d->qcn[k]=d->qcp[k]=0.; }
    for(int ax=0;ax<3;ax++){
        int n = ax==0?Nx:(ax==1?Ny:Nz);
        const double *c = ax==0?d->xs:(ax==1?d->ys:d->zs);
        size_t st = ax==0?1:(ax==1?(size_t)Nx:(size_t)Nx*Ny);
        int n1 = ax==0?Ny:Nx, n2 = ax==2?Ny:Nz;       /* the two other axes */
        #pragma omp parallel for collapse(2) schedule(static) if(N>PAR_MIN)
        for(int b=0;b<n2;b++) for(int a=0;a<n1;a++){
            size_t k0 = ax==0 ? IDX(0,a,b,Nx,Ny) : (ax==1 ? IDX(a,0,b,Nx,Ny) : IDX(a,b,0,Nx,Ny));
            for(int side=0;side<2;side++){            /* 0: wall on the low side, 1: high side */
                double w=0.; int have=0;
                for(int t=0;t<n;t++){
                    int q = side? n-1-t : t;
                    size_t k=k0+(size_t)q*st;
                    if(d->kind[k]!=K_SEMI){ have=0; continue; }
                    int qp = side? q+1 : q-1;         /* the neighbour we came from */
                    if(qp<0||qp>=n){
                        int g=d->gate_mask[k]; have=0;
                        if(g){ int f=d->cons[g-1].face; if(f/2==ax && (f%2)==side){ w=c[q]; have=1; } }
                    } else {
                        size_t kp=k0+(size_t)qp*st;
                        if(d->kind[kp]==K_INS){ w=0.5*(c[q]+c[qp]); have=1; }
                        else if(d->kind[kp]!=K_SEMI) have=0;
                    }
                    if(!have) continue;
                    double lo = q>0 ? 0.5*(c[q-1]+c[q]) : c[q];
                    double hi = q<n-1 ? 0.5*(c[q]+c[q+1]) : c[q];
                    double za = side? w-hi : lo-w, zb = side? w-lo : hi-w, zn=fabs(c[q]-w);
                    if(za<0.) za=0.;
                    if(zn<d->qd[k]) d->qd[k]=zn;
                    if(qc && za<WALL_REACH){
                        const Mat *m=MK(d,k);
                        d->qcn[k]+=log(fmax(mlda_factor(m,za,zb,0),1e-30));
                        d->qcp[k]+=log(fmax(mlda_factor(m,za,zb,1),1e-30));
                    }
                }
            }
        }
    }
    for(size_t k=0;k<N;k++){                      /* keep the exponents sane */
        if(d->qcn[k]<-60.) d->qcn[k]=-60.;
        if(d->qcp[k]<-60.) d->qcp[k]=-60.;
    }
}

/* Gate tunnelling paths on the tensor mesh.  From every semiconductor node
   (not a contact) with an insulator neighbour, march straight through the
   insulator; the path counts if it ends on the metal of a gate (bc 2
   Dirichlet node) - a semiconductor, another contact, the domain edge,
   more than TUN_LMAX of insulator or more than TUN_MAXL layers discard it.
   Layers: the semiconductor half cell (node -> interface midpoint), then
   one layer per run of equal insulator material, the last ending on the
   metal node (the electrostatic metal surface of the box method).  A Robin
   gate adds one path per gate node with the contact's stack.  Each path
   also gets its supply column: the semiconductor nodes inward from the
   interface over TUN_COL (stops at contacts and non-semiconductors). */
#define TUN_LMAX 30e-9
#define TUN_COL  5e-9
static int tun_col_push(Dev *d, size_t k, double l, double h){
    if(d->n_col>=d->cap_col){
        int cap = d->cap_col ? 2*d->cap_col : 8192;
        size_t *nk=realloc(d->col_k,cap*sizeof(size_t)); if(!nk) return 0; d->col_k=nk;
        double *nl=realloc(d->col_l,cap*sizeof(double)); if(!nl) return 0; d->col_l=nl;
        double *nh=realloc(d->col_h,cap*sizeof(double)); if(!nh) return 0; d->col_h=nh;
        double *nK=realloc(d->col_K,(size_t)cap*2*sizeof(double)); if(!nK) return 0; d->col_K=nK;
        d->cap_col=cap;
    }
    int q=d->n_col++;
    d->col_k[q]=k; d->col_l[q]=l; d->col_h[q]=h; d->col_K[2*q]=d->col_K[2*q+1]=0.;
    return 1;
}
static int tun_push(Dev *d, size_t k, double A, int c, int dir, int nl, const int *m, const double *t){
    if(d->n_tp>=d->cap_tp){
        int cap = d->cap_tp ? 2*d->cap_tp : 1024;
        size_t *nk=realloc(d->tp_k,cap*sizeof(size_t)); if(!nk) return 0; d->tp_k=nk;
        double *nA=realloc(d->tp_A,cap*sizeof(double)); if(!nA) return 0; d->tp_A=nA;
        int *nc=realloc(d->tp_c,cap*sizeof(int)); if(!nc) return 0; d->tp_c=nc;
        int *nd=realloc(d->tp_dir,cap*sizeof(int)); if(!nd) return 0; d->tp_dir=nd;
        int *nn=realloc(d->tp_nl,cap*sizeof(int)); if(!nn) return 0; d->tp_nl=nn;
        int *nm=realloc(d->tp_m,(size_t)cap*TUN_MAXL*sizeof(int)); if(!nm) return 0; d->tp_m=nm;
        double *nt=realloc(d->tp_t,(size_t)cap*TUN_MAXL*sizeof(double)); if(!nt) return 0; d->tp_t=nt;
        double *nG=realloc(d->tp_G,(size_t)cap*2*sizeof(double)); if(!nG) return 0; d->tp_G=nG;
        int *n0=realloc(d->tp_c0,cap*sizeof(int)); if(!n0) return 0; d->tp_c0=n0;
        int *n1=realloc(d->tp_cn,cap*sizeof(int)); if(!n1) return 0; d->tp_cn=n1;
        d->cap_tp=cap;
    }
    int p=d->n_tp;
    d->tp_k[p]=k; d->tp_A[p]=A; d->tp_c[p]=c; d->tp_dir[p]=dir; d->tp_nl[p]=nl;
    for(int i=0;i<TUN_MAXL;i++){ d->tp_m[p*TUN_MAXL+i]= i<nl?m[i]:0; d->tp_t[p*TUN_MAXL+i]= i<nl?t[i]:0.; }
    d->tp_G[2*p]=d->tp_G[2*p+1]=0.;
    /* supply column, inward (direction dir^1) */
    d->tp_c0[p]=d->n_col; int cnt=0;
    if(dir>=0){ int i,j,l; unidx(d,k,&i,&j,&l); size_t cur=k; double depth=0.; int in=dir^1, ax=dir/2;
      for(;;){
          int ii,jj,ll; unidx(d,cur,&ii,&jj,&ll);
          double L = ax==0? d->Lx[ii] : (ax==1? d->Ly[jj] : d->Lz[ll]);
          size_t kn; double h,ar; int more=nbr(d,ii,jj,ll,in,&kn,&h,&ar);
          if(!more || d->kind[kn]!=K_SEMI || d->con_mask[kn]) more=0;
          if(!tun_col_push(d,cur,L,more?h:0.)) return 0;
          cnt++; depth+=L;
          if(!more || depth>=TUN_COL) break;
          cur=kn;
      }
      (void)i; (void)j; (void)l; }
    d->tp_cn[p]=cnt;
    d->n_tp++;
    d->cons[c].ntp++;
    return 1;
}
static void tun_build(Dev *d){
    d->n_tp=0; d->n_col=0;
    for(int c=0;c<d->n_cons;c++) d->cons[c].ntp=0;
    if(!(d->models&32)) return;
    const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz;
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        if(d->kind[k]!=K_SEMI || d->con_mask[k]) continue;
        for(int dir=0;dir<6;dir++){
            size_t kb; double h,ar;
            if(!nbr(d,i,j,l,dir,&kb,&h,&ar) || d->kind[kb]!=K_INS || d->con_mask[kb]) continue;
            int m[TUN_MAXL]; double t[TUN_MAXL]; int nl=1, ok=0;
            m[0]=d->mat_idx[k]; t[0]=0.5*h;                /* semiconductor half cell */
            double tins=0., hin=h; size_t cur=kb; int ci,cj,cl; unidx(d,cur,&ci,&cj,&cl);
            for(int guard=0;guard<100000;guard++){
                size_t kn; double h2,ar2;
                if(!nbr(d,ci,cj,cl,dir,&kn,&h2,&ar2)) break;              /* domain edge */
                int mi=d->mat_idx[cur];
                int end = d->kind[kn]==K_METAL && d->con_mask[kn] && d->cons[d->con_mask[kn]-1].bc==2;
                double seg = 0.5*hin + (end ? h2 : 0.5*h2);  /* midpoint..midpoint (or ..metal) */
                if(nl>1 && m[nl-1]==mi) t[nl-1]+=seg;
                else { if(nl>=TUN_MAXL) break; m[nl]=mi; t[nl]=seg; nl++; }
                tins+=seg; if(tins>TUN_LMAX) break;
                if(end){ ok=d->con_mask[kn]; break; }
                if(d->kind[kn]!=K_INS || d->con_mask[kn]) break;
                cur=kn; hin=h2; unidx(d,cur,&ci,&cj,&cl);
            }
            if(ok && !tun_push(d,k,ar,ok-1,dir,nl,m,t)) return;
        }
    }
    for(int r=0;r<d->n_rob;r++){
        int c=d->rob_c[r]; const Con *cn=&d->cons[c];
        if(cn->snl<1) continue;
        if(!tun_push(d,d->rob_k[r],d->rob_A[r],c,cn->face,cn->snl,cn->smat,cn->st)) return;
    }
    if(d->verbose) fprintf(stderr,"semisim: %d tunnelling paths, %d column nodes\n",d->n_tp,d->n_col);
}

/* ---------------------------------------------- unstructured counterparts */
/* MLDA walls and interface distances supplied by the host: gwall holds per
   node two walls (distance of the node, range a..b of its cell; 1 = none) */
static void walls_graph(Dev *d, int qc){
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        d->qd[k]=1.; d->qcn[k]=d->qcp[k]=0.;
        if(d->kind[k]!=K_SEMI) continue;
        const Mat *m=MK(d,k);
        for(int w=0;w<2;w++){
            const double *g=d->gwall+6*k+3*w;
            if(!(g[0]<0.5)) continue;
            if(g[0]<d->qd[k]) d->qd[k]=g[0];
            if(qc && g[1]<WALL_REACH){
                d->qcn[k]+=log(fmax(mlda_factor(m,g[1],g[2],0),1e-30));
                d->qcp[k]+=log(fmax(mlda_factor(m,g[1],g[2],1),1e-30));
            }
        }
        if(d->qcn[k]<-60.) d->qcn[k]=-60.;
        if(d->qcp[k]<-60.) d->qcp[k]=-60.;
    }
}
static inline double g_proj(const Dev *d, size_t a, size_t b, const double *u){
    return (d->gpos[3*b]-d->gpos[3*a])*u[0]+(d->gpos[3*b+1]-d->gpos[3*a+1])*u[1]+(d->gpos[3*b+2]-d->gpos[3*a+2])*u[2];
}
/* the neighbour of cur best aligned with the unit vector u (cosine > cmin), not prev */
static long g_step(const Dev *d, size_t cur, long prev, const double *u, double cmin){
    long best=-1; double bc=cmin;
    for(int e=d->grp[cur];e<d->grp[cur+1];e++){
        size_t b=d->gcol[e]; if((long)b==prev) continue;
        double c=g_proj(d,cur,b,u)/d->gh[e];
        if(c>bc){ bc=c; best=(long)b; }
    }
    return best;
}
/* path + supply column on the graph: u = unit vector towards the gate, hw =
   distance node -> interface along u (0 for a face gate) */
static int tun_push_g(Dev *d, size_t k, double A, int c, const double *u, double hw,
                      int nl, const int *m, const double *t){
    if(!tun_push(d,k,A,c,-1,nl,m,t)) return 0;
    int p=d->n_tp-1;
    double v[3]={-u[0],-u[1],-u[2]};
    d->n_col=d->tp_c0[p];                      /* replace the (empty) tensor column */
    size_t cur=k; long prev=-1; double depth=0., hprev=2.*hw; int cnt=0;
    for(int guard=0;guard<10000;guard++){
        long nb=g_step(d,cur,prev,v,0.8);
        if(nb>=0 && (d->kind[nb]!=K_SEMI || d->con_mask[nb])) nb=-1;
        double h = nb>=0 ? g_proj(d,cur,(size_t)nb,v) : 0.;
        if(nb>=0 && !(h>0.)) nb=-1;
        double L = 0.5*hprev + (nb>=0 ? 0.5*h : 0.5*hprev);
        if(!(L>0.)) L=cbrt(d->gvol[cur]);
        if(!tun_col_push(d,cur,L,nb>=0?h:0.)) return 0;
        cnt++; depth+=L;
        if(nb<0 || depth>=TUN_COL) break;
        prev=(long)cur; cur=(size_t)nb; hprev=h;
    }
    d->tp_cn[p]=cnt;
    return 1;
}
static void tun_build_graph(Dev *d){
    d->n_tp=0; d->n_col=0;
    for(int c=0;c<d->n_cons;c++) d->cons[c].ntp=0;
    if(!(d->models&32)) return;
    size_t N=NN(d);
    for(size_t k=0;k<N;k++){
        if(d->kind[k]!=K_SEMI || d->con_mask[k]) continue;
        for(int e0=d->grp[k];e0<d->grp[k+1];e0++){
            size_t kb=d->gcol[e0];
            if(d->kind[kb]!=K_INS || d->con_mask[kb]) continue;
            double h=d->gh[e0], u[3];
            for(int q=0;q<3;q++) u[q]=(d->gpos[3*kb+q]-d->gpos[3*k+q])/h;
            int m[TUN_MAXL]; double t[TUN_MAXL]; int nl=1, ok=0;
            m[0]=d->mat_idx[k]; t[0]=0.5*h;
            double s0=0.5*h, scur=g_proj(d,k,kb,u), tins=0.;
            size_t cur=kb; long prev=(long)k;
            for(int guard=0;guard<100000;guard++){
                long nb=g_step(d,cur,prev,u,0.8);
                if(nb<0) break;
                double snb=g_proj(d,k,(size_t)nb,u);
                int end = d->kind[nb]==K_METAL && d->con_mask[nb] && d->cons[d->con_mask[nb]-1].bc==2;
                double s1 = end ? snb : 0.5*(scur+snb), seg=s1-s0;
                if(!(seg>0.)) break;
                int mi=d->mat_idx[cur];
                if(nl>1 && m[nl-1]==mi) t[nl-1]+=seg;
                else { if(nl>=TUN_MAXL) break; m[nl]=mi; t[nl]=seg; nl++; }
                tins+=seg; if(tins>TUN_LMAX) break;
                if(end){ ok=d->con_mask[nb]; break; }
                if(d->kind[nb]!=K_INS || d->con_mask[nb]) break;
                s0=s1; scur=snb; prev=(long)cur; cur=(size_t)nb;
            }
            if(ok && !tun_push_g(d,k,d->ga[e0],ok-1,u,0.5*h,nl,m,t)) return;
        }
    }
    for(int r=0;r<d->n_rob;r++){
        int c=d->rob_c[r]; const Con *cn=&d->cons[c];
        if(cn->snl<1) continue;
        double u[3]={0,0,0}; int f=cn->face; u[f/2] = (f&1) ? 1. : -1.;   /* outward normal of the face */
        if(!tun_push_g(d,d->rob_k[r],d->rob_A[r],c,u,0.,cn->snl,cn->smat,cn->st)) return;
    }
    if(d->verbose) fprintf(stderr,"semisim: %d tunnelling paths, %d column nodes (graph)\n",d->n_tp,d->n_col);
}

/* contact nodes: the host's node list (graph devices; shaped contacts on a
   tensor mesh) or the index ranges (tensor) */
#define FOR_CON_ANY(d,c,cn,...) do{ \
    if((d)->unstr || (d)->cnn[c]>0){ \
        for(int _q=0;_q<(d)->cnn[c];_q++){ size_t k=(size_t)(d)->cnodes[c][_q]; \
            double area=(d)->careas[c][_q]; (void)area; __VA_ARGS__ } \
    } else FOR_CON_NODES(d,cn,__VA_ARGS__); }while(0)

static int prepare(Dev *d){
    if(!d->allocated||d->n_mats<1) return 0;
    size_t N=NN(d);
    mat_prepare(d);
    for(size_t k=0;k<N;k++){
        int mi=d->mat_idx[k]; if(mi<0||mi>=d->n_mats) mi=d->mat_idx[k]=0;
        const Mat *m=&d->mats[mi];
        d->kind[k]= m->insulator ? K_INS : K_SEMI;
        d->eps[k] = m->eps_r*EPS0;
        d->bgn[k] = (!m->insulator && (d->models&1)) ? bgn_half(m,d->Nd[k]+d->Na[k]) : 0.;
        d->con_mask[k]=0; d->gate_mask[k]=0;
    }
    /* contacts: later contacts override earlier ones on shared nodes */
    int nrob=0;
    for(int c=0;c<d->n_cons;c++){
        Con *cn=&d->cons[c]; cn->A=0.;
        FOR_CON_ANY(d,c,cn,{
            if(cn->face<=5) cn->A+=area;
            if(cn->bc==2 && cn->face<=5 && d->kind[k]==K_SEMI){
                d->gate_mask[k]=c+1; nrob++;
            } else if(cn->bc==2 || d->kind[k]!=K_SEMI){
                d->kind[k]=K_METAL; d->con_mask[k]=c+1; d->gate_mask[k]=0;
            } else {
                d->con_mask[k]=c+1; d->gate_mask[k]=0;
            }
        });
    }
    free(d->rob_k); free(d->rob_A); free(d->rob_c);
    d->rob_k=NULL; d->rob_A=NULL; d->rob_c=NULL; d->n_rob=0;
    if(nrob>0){
        d->rob_k=malloc(nrob*sizeof(size_t)); d->rob_A=malloc(nrob*sizeof(double));
        d->rob_c=malloc(nrob*sizeof(int));
        if(!d->rob_k||!d->rob_A||!d->rob_c) return 0;
        for(int c=0;c<d->n_cons;c++){
            Con *cn=&d->cons[c];
            if(cn->bc!=2||cn->face>5) continue;
            FOR_CON_ANY(d,c,cn,{
                if(d->gate_mask[k]==c+1 && !d->con_mask[k] && d->n_rob<nrob){
                    d->rob_k[d->n_rob]=k; d->rob_A[d->n_rob]=area; d->rob_c[d->n_rob]=c;
                    d->n_rob++;
                }
            });
        }
    }
    /* box contacts: area = interface with non-contact semiconductor */
    for(int c=0;c<d->n_cons;c++){
        Con *cn=&d->cons[c]; if(cn->face<=5) continue;
        double A=0.;
        if(d->unstr){
            for(int q=0;q<d->cnn[c];q++){
                size_t k=(size_t)d->cnodes[c][q]; if(d->con_mask[k]!=c+1) continue;
                for(int e=d->grp[k];e<d->grp[k+1];e++){
                    size_t kb=d->gcol[e];
                    if(d->kind[kb]==K_SEMI && d->con_mask[kb]!=c+1) A+=d->ga[e];
                }
            }
        } else if(d->cnn[c]>0){
            for(int q=0;q<d->cnn[c];q++){
                size_t k=(size_t)d->cnodes[c][q]; if(d->con_mask[k]!=c+1) continue;
                int i,j,l; unidx(d,k,&i,&j,&l);
                for(int dir=0;dir<6;dir++){
                    size_t kb; double h,ar;
                    if(nbr(d,i,j,l,dir,&kb,&h,&ar) && d->kind[kb]==K_SEMI && d->con_mask[kb]!=c+1) A+=ar;
                }
            }
        } else FOR_CON_NODES(d,cn,{
            if(d->con_mask[k]==c+1) for(int dir=0;dir<6;dir++){
                size_t kb; double h,ar;
                if(nbr(d,i,j,l,dir,&kb,&h,&ar) && d->kind[kb]==K_SEMI && d->con_mask[kb]!=c+1) A+=ar;
            }
        });
        cn->A=A;
    }
    if(d->unstr){
        walls_graph(d,(d->models&8)!=0);
        g_lsq_setup(d);
        tun_build_graph(d);
    } else {
        walls_tensor(d,(d->models&8)!=0);
        tun_build(d);
    }
    d->dirty=0;
    return 1;
}

/* Dirichlet values for contact nodes at the currently applied voltages (Vapp) */
static void apply_bc(Dev *d, const double *Vapp){
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        int c=d->con_mask[k]; if(!c) continue;
        const Con *cn=&d->cons[c-1]; double V=Vapp[c-1];
        const Mat *m=&d->mats[d->mat_idx[k]];
        if(d->kind[k]==K_METAL){
            double pm = cn->phi_m>0. ? cn->phi_m : d->Cref;
            d->phi[k]=V+d->Cref-pm; d->phi_n[k]=d->phi_p[k]=V; d->n[k]=d->p[k]=0.;
            continue;
        }
        if(cn->bc==1){
            double pm = cn->phi_m>0. ? cn->phi_m : d->Cref;
            d->phi[k]=V+d->Cref-pm;
        } else {
            double D=d->Nd[k]-d->Na[k];
            d->phi[k]=V+m->VT*(asinh(D/(2.*niK(d,k)))-0.5*(lncK(d,k)-lnvK(d,k)));
        }
        d->phi_n[k]=d->phi_p[k]=V;
        d->n[k]=nK(d,k,d->phi[k],V); d->p[k]=pK(d,k,d->phi[k],V);
    }
}

/* Local charge-neutral equilibrium start (all quasi-Fermi potentials 0) */
static void init_eq(Dev *d){
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        const Mat *m=&d->mats[d->mat_idx[k]];
        d->phi_n[k]=d->phi_p[k]=0.;
        if(d->kind[k]==K_SEMI){
            double D=d->Nd[k]-d->Na[k];
            d->phi[k]=m->VT*(asinh(D/(2.*niK(d,k)))-0.5*(lncK(d,k)-lnvK(d,k)));
            d->n[k]=nK(d,k,d->phi[k],0.); d->p[k]=pK(d,k,d->phi[k],0.);
            d->mu_n[k]=mob_low(m,d->Nd[k]+d->Na[k],0);
            d->mu_p[k]=mob_low(m,d->Nd[k]+d->Na[k],1);
        } else { d->phi[k]=0.; d->n[k]=d->p[k]=0.; d->mu_n[k]=d->mu_p[k]=0.; }
        d->R_tot[k]=0.;
    }
}

/* ------------------------------------------------------ linear algebra
   7-point stencil matrices on tensor grids.
   Smoother / fine preconditioner: chunked D-ILU(0) (exactly ILU(0) for a
   7-point stencil in natural ordering; IC(0) when A is symmetric). Chunks
   are contiguous index ranges, one per thread (block Jacobi across chunks).
   Preconditioner: multigrid V(1,1) cycle with 2x2x2 aggregation. For a
   7-point stencil an aggregation coarse operator R A P is again a 7-point
   stencil, so every level reuses the same kernels; the coarsest level
   (<= 64 unknowns) is solved by dense LU.
   Poisson (SPD): P = piecewise constant, R = P^T -> symmetric V-cycle (CG).
   Continuity (conservation form, unknown u = n/s): R sums the UNSCALED flux
   residuals (left null vector 1) and P is Slotboom-weighted,
   P_kI = exp((phi_k - phi_I)/VT) with phi = phi_n (electrons) or -phi_p
   (holes) and phi_I the aggregate maximum: a coarse correction shifts the
   quasi-Fermi level of an aggregate uniformly, which is the near-null space
   of the drift-diffusion operator even when the current iterate is far from
   it. At convergence inside an aggregate this reduces to constant P. */
static double MG_ALPHA=1.0;

static int n_chunks_N(const Dev *d, size_t N){
    if(N<(size_t)PAR_MIN||d->nthreads<2) return 1;
    int nc=d->nthreads; if((size_t)nc*4096>N) nc=(int)(N/4096); if(nc<1)nc=1;
    return nc;
}

static void sten_free(Sten *A){
    free(A->c);free(A->xm);free(A->xp);free(A->ym);free(A->yp);free(A->zm);free(A->zp);
    memset(A,0,sizeof(*A));
}
static int sten_alloc(Sten *A, size_t N){
    A->c=calloc(N,8);A->xm=calloc(N,8);A->xp=calloc(N,8);A->ym=calloc(N,8);
    A->yp=calloc(N,8);A->zm=calloc(N,8);A->zp=calloc(N,8);
    return A->c&&A->xm&&A->xp&&A->ym&&A->yp&&A->zm&&A->zp;
}

static void mg_free(Dev *d){
    for(int L=0;L<MAXLEV;L++){
        Lev *v=&d->lev[L];
        if(v->own){ sten_free(&v->A); free(v->ilu); free(v->x); free(v->b); free(v->ph);
                    free(v->rp); free(v->cj); free(v->rv); free(v->va); free(v->dg); }
        free(v->r); free(v->pw); free(v->LU); free(v->piv);
        free(v->agg); free(v->mrp); free(v->mem); free(v->emap); free(v->wg);
        memset(v,0,sizeof(*v));
    }
    d->nlev=0;
}

static int mg_alloc(Dev *d){
    mg_free(d);
    Lev *f=&d->lev[0];
    f->Nx=d->Nx; f->Ny=d->Ny; f->Nz=d->Nz; f->N=NN(d);
    f->A=d->A; f->ilu=d->ilu; f->r=calloc(f->N,8); f->pw=calloc(f->N,8);
    if(!f->r||!f->pw) return 0;
    int nl=1;
    while(nl<MAXLEV){
        f=&d->lev[nl-1];
        if(f->N<=64) break;
        int cx=f->Nx>2, cy=f->Ny>2, cz=f->Nz>2;
        if(!cx&&!cy&&!cz) break;
        f->cx=cx; f->cy=cy; f->cz=cz;
        Lev *c=&d->lev[nl];
        c->Nx=cx?(f->Nx+1)/2:f->Nx; c->Ny=cy?(f->Ny+1)/2:f->Ny; c->Nz=cz?(f->Nz+1)/2:f->Nz;
        c->N=(size_t)c->Nx*c->Ny*c->Nz; c->own=1;
        if(!sten_alloc(&c->A,c->N)) return 0;
        c->ilu=calloc(c->N,8); c->x=calloc(c->N,8); c->b=calloc(c->N,8); c->r=calloc(c->N,8);
        c->pw=calloc(c->N,8); c->ph=calloc(c->N,8);
        if(!c->ilu||!c->x||!c->b||!c->r||!c->pw||!c->ph) return 0;
        nl++;
    }
    Lev *z=&d->lev[nl-1];
    z->LU=calloc(z->N*z->N,8); z->piv=calloc(z->N,sizeof(int));
    if(!z->LU||!z->piv) return 0;
    d->nlev=nl;
    return 1;
}

static void ilu_factor_L(const Dev *d, Lev *v){
    const int Nx=v->Nx,Ny=v->Ny; const size_t NxNy=(size_t)Nx*Ny, N=v->N;
    const Sten *A=&v->A; double *D=v->ilu;
    int nch=n_chunks_N(d,N);
    #pragma omp parallel for schedule(static,1) if(nch>1)
    for(int c=0;c<nch;c++){
        size_t k0=N*(size_t)c/nch, k1=N*(size_t)(c+1)/nch;
        int i=(int)(k0%Nx), j=(int)((k0/Nx)%Ny);
        for(size_t k=k0;k<k1;k++){
            double dk=A->c[k];
            if(i>0 && k>k0)                 dk-=A->xm[k]*A->xp[k-1]   /D[k-1];
            if(j>0 && k>=k0+(size_t)Nx)     dk-=A->ym[k]*A->yp[k-Nx]  /D[k-Nx];
            if(k>=k0+NxNy && k>=NxNy)       dk-=A->zm[k]*A->zp[k-NxNy]/D[k-NxNy];
            if(!(dk>1e-8*fabs(A->c[k]))) dk=A->c[k];
            if(!(dk!=0.)) dk=1.;
            D[k]=dk;
            if(++i==Nx){ i=0; if(++j==Ny) j=0; }
        }
    }
}

/* z = M^-1 r ; in place allowed (z == r) */
static void ilu_apply_L(const Dev *d, const Lev *v, const double *r, double *z){
    const int Nx=v->Nx,Ny=v->Ny; const size_t NxNy=(size_t)Nx*Ny, N=v->N;
    const Sten *A=&v->A; const double *D=v->ilu;
    int nch=n_chunks_N(d,N);
    #pragma omp parallel for schedule(static,1) if(nch>1)
    for(int c=0;c<nch;c++){
        size_t k0=N*(size_t)c/nch, k1=N*(size_t)(c+1)/nch;
        int i=(int)(k0%Nx), j=(int)((k0/Nx)%Ny);
        for(size_t k=k0;k<k1;k++){
            double s=r[k];
            if(i>0 && k>k0)             s-=A->xm[k]*z[k-1];
            if(j>0 && k>=k0+(size_t)Nx) s-=A->ym[k]*z[k-Nx];
            if(k>=k0+NxNy)              s-=A->zm[k]*z[k-NxNy];
            z[k]=s/D[k];
            if(++i==Nx){ i=0; if(++j==Ny) j=0; }
        }
        size_t kl=k1-1;
        i=(int)(kl%Nx); j=(int)((kl/Nx)%Ny);
        for(size_t kk=k1;kk>k0;kk--){
            size_t k=kk-1; double s=0.;
            if(i<Nx-1 && k+1<k1)            s+=A->xp[k]*z[k+1];
            if(j<Ny-1 && k+(size_t)Nx<k1)   s+=A->yp[k]*z[k+Nx];
            if(k+NxNy<k1)                   s+=A->zp[k]*z[k+NxNy];
            z[k]-=s/D[k];
            if(--i<0){ i=Nx-1; if(--j<0) j=Ny-1; }
        }
    }
}

/* y = A x  (res: y = b - A x when b != NULL) */
static void matvec_L(const Dev *d, const Lev *v, const double *x, const double *b, double *y){
    const int Nx=v->Nx,Ny=v->Ny,Nz=v->Nz; const size_t NxNy=(size_t)Nx*Ny;
    const Sten *A=&v->A;
    #pragma omp parallel for collapse(2) schedule(static) if(v->N>PAR_MIN)
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++){
        size_t k=IDX(0,j,l,Nx,Ny);
        for(int i=0;i<Nx;i++,k++){
            double s=A->c[k]*x[k];
            if(i>0)    s+=A->xm[k]*x[k-1];
            if(i<Nx-1) s+=A->xp[k]*x[k+1];
            if(j>0)    s+=A->ym[k]*x[k-Nx];
            if(j<Ny-1) s+=A->yp[k]*x[k+Nx];
            if(l>0)    s+=A->zm[k]*x[k-NxNy];
            if(l<Nz-1) s+=A->zp[k]*x[k+NxNy];
            y[k]= b ? b[k]-s : s;
        }
    }
    (void)d;
}

#define AGG_RANGE(f,I,J,L) \
    int i0=(f)->cx?2*(I):(I), i1=(f)->cx?2*(I)+1:(I); if(i1>(f)->Nx-1) i1=(f)->Nx-1; \
    int j0=(f)->cy?2*(J):(J), j1=(f)->cy?2*(J)+1:(J); if(j1>(f)->Ny-1) j1=(f)->Ny-1; \
    int l0=(f)->cz?2*(L):(L), l1=(f)->cz?2*(L)+1:(L); if(l1>(f)->Nz-1) l1=(f)->Nz-1;

/* Coarse operator P^T W A P (W = diag weights w, or identity). With the
   continuity rows scaled by 1/D_a, w = D restores the conservation form, so
   the coarse operator sums fluxes (Petrov-Galerkin, left null vector = 1). */
static void galerkin(const Lev *f, Lev *c, const double *w){
    const Sten *A=&f->A;
    #pragma omp parallel for collapse(2) schedule(static) if(c->N>PAR_MIN/8)
    for(int L=0;L<c->Nz;L++) for(int J=0;J<c->Ny;J++) for(int I=0;I<c->Nx;I++){
        AGG_RANGE(f,I,J,L)
        double cc=0,xm=0,xp=0,ym=0,yp=0,zm=0,zp=0;
        for(int l=l0;l<=l1;l++) for(int j=j0;j<=j1;j++) for(int i=i0;i<=i1;i++){
            size_t k=IDX(i,j,l,f->Nx,f->Ny); double a, wk= w? w[k] : 1.;
            const double *pw=f->pw; const size_t NxNy=(size_t)f->Nx*f->Ny;
            cc+=wk*A->c[k]*pw[k];
            if(i>0)      { a=wk*A->xm[k]*pw[k-1];    if(i-1>=i0) cc+=a; else xm+=a; }
            if(i<f->Nx-1){ a=wk*A->xp[k]*pw[k+1];    if(i+1<=i1) cc+=a; else xp+=a; }
            if(j>0)      { a=wk*A->ym[k]*pw[k-f->Nx]; if(j-1>=j0) cc+=a; else ym+=a; }
            if(j<f->Ny-1){ a=wk*A->yp[k]*pw[k+f->Nx]; if(j+1<=j1) cc+=a; else yp+=a; }
            if(l>0)      { a=wk*A->zm[k]*pw[k-NxNy];  if(l-1>=l0) cc+=a; else zm+=a; }
            if(l<f->Nz-1){ a=wk*A->zp[k]*pw[k+NxNy];  if(l+1<=l1) cc+=a; else zp+=a; }
        }
        size_t K=IDX(I,J,L,c->Nx,c->Ny);
        /* an aggregate made only of identity rows (contacts, insulators)
           carries no equation: keep it as an identity row */
        if(cc==0.&&xm==0.&&xp==0.&&ym==0.&&yp==0.&&zm==0.&&zp==0.) cc=1.;
        c->A.c[K]=cc; c->A.xm[K]=xm; c->A.xp[K]=xp; c->A.ym[K]=ym;
        c->A.yp[K]=yp; c->A.zm[K]=zm; c->A.zp[K]=zp;
    }
}
static void restrict_sum(const Lev *f, const Lev *c, const double *r, double *bc, const double *w){
    #pragma omp parallel for collapse(2) schedule(static) if(c->N>PAR_MIN/8)
    for(int L=0;L<c->Nz;L++) for(int J=0;J<c->Ny;J++) for(int I=0;I<c->Nx;I++){
        AGG_RANGE(f,I,J,L)
        double s=0.;
        for(int l=l0;l<=l1;l++) for(int j=j0;j<=j1;j++) for(int i=i0;i<=i1;i++){
            size_t k=IDX(i,j,l,f->Nx,f->Ny); s+= w ? w[k]*r[k] : r[k]; }
        bc[IDX(I,J,L,c->Nx,c->Ny)]=s;
    }
}
static void prolong_add(const Lev *f, const Lev *c, const double *xc, double *x){
    #pragma omp parallel for collapse(2) schedule(static) if(f->N>PAR_MIN)
    for(int l=0;l<f->Nz;l++) for(int j=0;j<f->Ny;j++){
        int L=f->cz?l>>1:l, J=f->cy?j>>1:j;
        size_t k=IDX(0,j,l,f->Nx,f->Ny);
        for(int i=0;i<f->Nx;i++,k++){
            int I=f->cx?i>>1:i;
            x[k]+=MG_ALPHA*f->pw[k]*xc[IDX(I,J,L,c->Nx,c->Ny)];
        }
    }
}

static void dense_factor(Lev *v){
    int n=(int)v->N; double *M=v->LU; const Sten *A=&v->A;
    memset(M,0,(size_t)n*n*8);
    for(int l=0;l<v->Nz;l++) for(int j=0;j<v->Ny;j++) for(int i=0;i<v->Nx;i++){
        int k=(int)IDX(i,j,l,v->Nx,v->Ny), NxNy=v->Nx*v->Ny;
        M[k*n+k]=A->c[k];
        if(i>0)       M[k*n+k-1]   =A->xm[k];
        if(i<v->Nx-1) M[k*n+k+1]   =A->xp[k];
        if(j>0)       M[k*n+k-v->Nx]=A->ym[k];
        if(j<v->Ny-1) M[k*n+k+v->Nx]=A->yp[k];
        if(l>0)       M[k*n+k-NxNy]=A->zm[k];
        if(l<v->Nz-1) M[k*n+k+NxNy]=A->zp[k];
    }
    for(int c=0;c<n;c++){
        int p=c; double mx=fabs(M[c*n+c]);
        for(int r=c+1;r<n;r++) if(fabs(M[r*n+c])>mx){ mx=fabs(M[r*n+c]); p=r; }
        v->piv[c]=p;
        if(p!=c) for(int q=0;q<n;q++){ double t=M[c*n+q]; M[c*n+q]=M[p*n+q]; M[p*n+q]=t; }
        double dg=M[c*n+c]; if(!(fabs(dg)>1e-300)) dg=M[c*n+c]=1e-300;
        for(int r=c+1;r<n;r++){
            double f=M[r*n+c]/dg; M[r*n+c]=f;
            if(f!=0.) for(int q=c+1;q<n;q++) M[r*n+q]-=f*M[c*n+q];
        }
    }
}
static void dense_solve(const Lev *v, const double *b, double *x){
    int n=(int)v->N; const double *M=v->LU;
    for(int k=0;k<n;k++) x[k]=b[k];
    for(int c=0;c<n;c++){ int p=v->piv[c]; if(p!=c){ double t=x[c]; x[c]=x[p]; x[p]=t; } }
    for(int r=0;r<n;r++){ double s=x[r]; for(int q=0;q<r;q++) s-=M[r*n+q]*x[q]; x[r]=s; }
    for(int r=n-1;r>=0;r--){ double s=x[r]; for(int q=r+1;q<n;q++) s-=M[r*n+q]*x[q]; x[r]=s/M[r*n+r]; }
}

#define PH_NONE (-1e300)
/* level-0 reference potential (units of VT) of node k */
static inline double ph0(const Dev *d, size_t k){
    if(d->mg_qsign==0) return 0.;
    if(d->kind[k]!=K_SEMI || d->con_mask[k]) return PH_NONE;
    double VT=d->mats[d->ref_mat].VT;
    return d->mg_qsign>0 ? d->phi_n[k]/VT : -d->phi_p[k]/VT;
}
/* Slotboom prolongation weights for every level (all 1 for mg_qsign==0) */
static void mg_weights(Dev *d){
    for(int L=0;L<d->nlev-1;L++){
        Lev *f=&d->lev[L], *c=&d->lev[L+1];
        if(d->mg_qsign==0){
            for(size_t k=0;k<f->N;k++) f->pw[k]=1.;
            continue;
        }
        #pragma omp parallel for collapse(2) schedule(static) if(c->N>PAR_MIN/8)
        for(int LL=0;LL<c->Nz;LL++) for(int J=0;J<c->Ny;J++) for(int I=0;I<c->Nx;I++){
            AGG_RANGE(f,I,J,LL)
            double mx=PH_NONE;
            for(int l=l0;l<=l1;l++) for(int j=j0;j<=j1;j++) for(int i=i0;i<=i1;i++){
                size_t k=IDX(i,j,l,f->Nx,f->Ny);
                double p= L==0 ? ph0(d,k) : f->ph[k];
                if(p>mx) mx=p;
            }
            c->ph[IDX(I,J,LL,c->Nx,c->Ny)]=mx;
            for(int l=l0;l<=l1;l++) for(int j=j0;j<=j1;j++) for(int i=i0;i<=i1;i++){
                size_t k=IDX(i,j,l,f->Nx,f->Ny);
                double p= L==0 ? ph0(d,k) : f->ph[k];
                double e= (p<=0.5*PH_NONE) ? -800. : p-mx;
                f->pw[k]= e<-700. ? 0. : exp(e);
            }
        }
    }
}

/* graph-mode counterparts (part2c) */
static void gmg_setup(Dev *d);
static void gvcycle(Dev *d, int L, const double *b, double *x);
static void gmatvec_L(const Dev *d, const Lev *v, const double *x, const double *b, double *y);
static double gfloor(const Dev *d, const double *b, const double *x);

/* Build the hierarchy from the current level-0 matrix (d->A) */
static void mg_setup(Dev *d){
    if(d->unstr){ gmg_setup(d); return; }
    double t0=wtime();
    mg_weights(d);
    for(int L=0;L<d->nlev-1;L++){
        galerkin(&d->lev[L],&d->lev[L+1], L==0 ? d->mg_w : NULL);
        ilu_factor_L(d,&d->lev[L]);
    }
    dense_factor(&d->lev[d->nlev-1]);
    d->t_mg+=wtime()-t0;
    if(d->verbose>5){
        for(int L=0;L<d->nlev;L++){
            Lev *v=&d->lev[L]; double mn=1e300,mx=0; long zr=0,nan=0;
            for(size_t k=0;k<v->N;k++){ double c=v->A.c[k];
                if(c!=c) nan++; if(fabs(c)<mn) mn=fabs(c); if(fabs(c)>mx) mx=fabs(c);
                if(c==0.&&v->A.xm[k]==0.&&v->A.xp[k]==0.) zr++; }
            fprintf(stderr,"        lev%d %dx%dx%d diag[%.2e,%.2e] zero-rows=%ld nan=%ld\n",L,v->Nx,v->Ny,v->Nz,mn,mx,zr,nan);
        }
        Lev *z=&d->lev[d->nlev-1]; double pm=1e300; int n=(int)z->N;
        for(int c=0;c<n;c++) if(fabs(z->LU[c*n+c])<pm) pm=fabs(z->LU[c*n+c]);
        fprintf(stderr,"        coarsest min pivot %.3e\n",pm);
    }
}

static int mg_off=0;
static void vcycle(Dev *d, int L, const double *b, double *x){
    Lev *v=&d->lev[L];
    if(mg_off && L==0 && d->nlev>1){ ilu_apply_L(d,v,b,x); return; }
    if(L==d->nlev-1){ dense_solve(v,b,x); return; }
    Lev *c=&d->lev[L+1]; size_t N=v->N;
    ilu_apply_L(d,v,b,x);                       /* pre-smoothing from x=0 */
    matvec_L(d,v,x,b,v->r);                     /* r = b - A x           */
    restrict_sum(v,c,v->r,c->b, L==0 ? d->mg_w : NULL);
    vcycle(d,L+1,c->b,c->x);
    prolong_add(v,c,c->x,x);
    matvec_L(d,v,x,b,v->r);                     /* post-smoothing        */
    ilu_apply_L(d,v,v->r,v->r);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++) x[k]+=v->r[k];
}

static inline void precond(Dev *d, const double *r, double *z){
    if(d->unstr) gvcycle(d,0,r,z); else vcycle(d,0,r,z);
}
static inline void matvec(Dev *d, const double *x, double *y){
    if(d->unstr) gmatvec_L(d,&d->lev[0],x,NULL,y); else matvec_L(d,&d->lev[0],x,NULL,y);
}

static double dot(const Dev *d, const double *a, const double *b){
    size_t N=NN(d); double s=0.;
    #pragma omp parallel for reduction(+:s) schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++) s+=a[k]*b[k];
    return s;
}

/* Preconditioned CG, x in/out. Stops on the M^-1 norm of the residual. */
static int pcg(Dev *d, const double *b, double *x, int maxit, double rtol){
    size_t N=NN(d);
    double *r=d->wk[0],*z=d->wk[1],*p=d->wk[2],*q=d->wk[3];
    mg_setup(d);
    matvec(d,x,q);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++) r[k]=b[k]-q[k];
    precond(d,r,z);
    double rz=dot(d,r,z), rz0=rz;
    if(!(rz0>0.)) return 0;
    memcpy(p,z,N*sizeof(double));
    int it;
    for(it=1;it<=maxit;it++){
        matvec(d,p,q);
        double pq=dot(d,p,q); if(!(pq>0.)){ if(d->verbose>3) fprintf(stderr,"      pcg breakdown pq=%.3e it=%d\n",pq,it); break; }
        double al=rz/pq;
        #pragma omp parallel for schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){ x[k]+=al*p[k]; r[k]-=al*q[k]; }
        precond(d,r,z);
        double rzn=dot(d,r,z);
        if(!(rzn>rtol*rtol*rz0)) break;
        double be=rzn/rz; rz=rzn;
        #pragma omp parallel for schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++) p[k]=z[k]+be*p[k];
    }
    if(d->verbose>3) fprintf(stderr,"      pcg its=%d rz0=%.3e rz=%.3e\n",it,rz0,rz);
    d->pcg_its+=it;
    return it;
}

/* Right-preconditioned BiCGSTAB, x in/out.
   Stops at ||r|| <= max(rtol*||r0||, atol*sqrt(N)); returns final ||r||. */
static double bicgstab_core(Dev *d, const double *b, double *x, int maxit,
                            double rtol, double atol, int *its, double *r0out);
static double bicgstab(Dev *d, const double *b, double *x, int maxit,
                       double rtol, double atol, int *its){
    size_t N=NN(d);
    double *x0=d->xsave;
    memcpy(x0,x,N*sizeof(double));
    double r0, rn=bicgstab_core(d,b,x,maxit,rtol,atol,its,&r0);
    if(!(rn<=fmax(r0,1e-300)) && d->nlev>1){   /* diverged / NaN: ILU only */
        if(d->verbose) fprintf(stderr,"semisim: MG-BiCGSTAB diverged (%.2e > %.2e), retry ILU\n",rn,r0);
        memcpy(x,x0,N*sizeof(double));
        int its2; mg_off=1;
        rn=bicgstab_core(d,b,x,4*maxit,rtol,atol,&its2,&r0);
        mg_off=0; *its+=its2;
    }
    /* still worse than the start (or not finite): keep the start, report
       failure - a diverged iterate must never enter the Gummel state */
    if(!(rn<=fmax(r0,1e-300)) || !isfinite(rn)){ memcpy(x,x0,N*sizeof(double)); return INFINITY; }
    return rn;
}
static double bicgstab_core(Dev *d, const double *b, double *x, int maxit,
                            double rtol, double atol, int *its, double *r0out){
    size_t N=NN(d);
    double *r=d->wk[0],*rh=d->wk[1],*p=d->wk[2],*v=d->wk[3],
           *ph=d->wk[4],*sh=d->wk[5],*t=d->wk[6];
    mg_setup(d);
    matvec(d,x,t);
    double rr=0.;
    #pragma omp parallel for reduction(+:rr) schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){ r[k]=b[k]-t[k]; rh[k]=r[k]; p[k]=0.; v[k]=0.; rr+=r[k]*r[k]; }
    /* rounding floor of this system: every row term carries a relative
       error ~atol (the scaled entries come from exp() of arguments up to
       ~100), so residuals below atol*||sum_b |A_ab x_b| + |b_a| || are noise
       and must not drive the solve (floating regions would amplify them). */
    double fl=0.;
    if(d->unstr) fl=gfloor(d,b,x);
    else {
        const Lev *v0=&d->lev[0]; const Sten *A=&v0->A;
        const int Nx=v0->Nx,Ny=v0->Ny,Nz=v0->Nz; const size_t NxNy=(size_t)Nx*Ny;
        #pragma omp parallel for collapse(2) reduction(+:fl) schedule(static) if(N>PAR_MIN)
        for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++){
            size_t k=IDX(0,j,l,Nx,Ny);
            for(int i=0;i<Nx;i++,k++){
                double s=fabs(A->c[k]*x[k])+fabs(b[k]);
                if(i>0)    s+=fabs(A->xm[k]*x[k-1]);
                if(i<Nx-1) s+=fabs(A->xp[k]*x[k+1]);
                if(j>0)    s+=fabs(A->ym[k]*x[k-Nx]);
                if(j<Ny-1) s+=fabs(A->yp[k]*x[k+Nx]);
                if(l>0)    s+=fabs(A->zm[k]*x[k-NxNy]);
                if(l<Nz-1) s+=fabs(A->zp[k]*x[k+NxNy]);
                fl+=s*s;
            }
        }
    }
    double r0=sqrt(rr), tgt=fmax(rtol*r0, atol*sqrt(fl));
    if(d->verbose>3) fprintf(stderr,"      r0=%.3e tgt=%.3e\n",r0,tgt);
    *its=0; *r0out=r0;
    if(d->verbose>6){
        size_t km=0; double mx=0; for(size_t k=0;k<N;k++) if(fabs(r[k])>mx){ mx=fabs(r[k]); km=k; }
        int i,j,l; unidx(d,km,&i,&j,&l);
        fprintf(stderr,"      max|r|=%.3e at (%d,%d,%d) s=%.3e n=%.3e p=%.3e phi=%.4f phin=%.4f phip=%.4f D=%.3e\n",
            mx,i,j,l,d->sc[km],d->n[km],d->p[km],d->phi[km],d->phi_n[km],d->phi_p[km],d->dsc[km]);
        for(int q=-3;q<=3;q++){ long kk=(long)km+q; if(kk<0||kk>=(long)N) continue;
            fprintf(stderr,"        k%+d: r=%.3e s=%.3e phi=%.4f phin=%.4f phip=%.4f\n",q,r[kk],d->sc[kk],d->phi[kk],d->phi_n[kk],d->phi_p[kk]); }
    }
    if(!(r0>tgt)) return r0;
    double rho=1.,al=1.,om=1., rbest=r0; int stall=0, it;
    for(it=1;it<=maxit;it++){
        double rho1=dot(d,rh,r);
        if(fabs(rho1)<1e-300){               /* breakdown: restart shadow */
            #pragma omp parallel for schedule(static) if(N>PAR_MIN)
            for(size_t k=0;k<N;k++){ rh[k]=r[k]; p[k]=0.; v[k]=0.; }
            rho1=dot(d,r,r); rho=al=om=1.;
            if(!(rho1>0.)) break;
        }
        double be=(rho1/rho)*(al/om);
        #pragma omp parallel for schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++) p[k]=r[k]+be*(p[k]-om*v[k]);
        precond(d,p,ph); matvec(d,ph,v);
        double rv=dot(d,rh,v); if(fabs(rv)<1e-300) break;
        al=rho1/rv;
        double ss=0.;
        #pragma omp parallel for reduction(+:ss) schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){ r[k]-=al*v[k]; ss+=r[k]*r[k]; }
        if(sqrt(ss)<=tgt){
            #pragma omp parallel for schedule(static) if(N>PAR_MIN)
            for(size_t k=0;k<N;k++) x[k]+=al*ph[k];
            rr=ss; break;
        }
        precond(d,r,sh); matvec(d,sh,t);
        double tt=0.,ts=0.;
        #pragma omp parallel for reduction(+:tt,ts) schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){ tt+=t[k]*t[k]; ts+=t[k]*r[k]; }
        om = tt>0. ? ts/tt : 0.;
        rr=0.;
        #pragma omp parallel for reduction(+:rr) schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){ x[k]+=al*ph[k]+om*sh[k]; r[k]-=om*t[k]; rr+=r[k]*r[k]; }
        double rn=sqrt(rr);
        if(d->verbose>5) fprintf(stderr,"        it %d rn=%.3e om=%.3e al=%.3e\n",it,rn,om,al);
        if(rn<=tgt || om==0.) break;
        if(rn<0.99*rbest){ rbest=rn; stall=0; } else if(++stall>40) break;
        rho=rho1;
    }
    if(it>maxit) it=maxit;
    d->bicg_its+=it; *its=it;
    return sqrt(rr);
}

/* ================================================ unstructured (graph) mode
   Linear algebra on a general symmetric-pattern sparse matrix (CSR: row
   pointer rp, column cj, reverse entry rv, off-diagonal values va, diagonal
   dg).  Used when the device comes from a triangular-prism (or any other
   box-method) mesh.  Same scheme as the tensor path:
     smoother   chunked D-ILU(0) (block Jacobi across contiguous chunks),
     multigrid  greedy aggregation (a root and its free neighbours), piece-
                wise constant (Poisson) or Slotboom-weighted (continuity)
                prolongation, Galerkin P^T W A P, dense LU at the coarsest.
   The aggregation depends only on the graph and is built once per mesh. */
#define GMG_COARSE 160           /* dense LU at or below this many unknowns */

static int g_chunks(const Dev *d, size_t N){
    if(N<(size_t)PAR_MIN||d->nthreads<2) return 1;
    int nc=d->nthreads; if((size_t)nc*4096>N) nc=(int)(N/4096); if(nc<1)nc=1;
    return nc;
}

/* Aggregation of level f by strength-based pairwise matching (as in Notay's
   aggregation AMG): every node joins its strongest free neighbour - the
   geometric coupling w = dual area / length, summed over merged nodes - if
   that coupling is at least GMG_BETA of its strongest one, else it stays
   alone.  GMG_PASSES rounds of matching give aggregates of up to 2^passes
   nodes that follow the strong direction of anisotropic cells (inversion
   layers, thin prism layers) instead of the plain neighbourhood, which on
   such meshes made coarse corrections weak and left smooth errors in the
   continuity solves.  agg[i] = aggregate, returns their number. */
static int GMG_PASSES=3; static double GMG_BETA=0.5;    /* SEMISIM_GPASS / SEMISIM_GBETA override */
static int g_pair(int N, const int *rp, const int *cj, const double *w, int *agg){
    for(int i=0;i<N;i++) agg[i]=-1;
    int na=0;
    for(int i=0;i<N;i++){
        if(agg[i]>=0) continue;
        double mx=0.;
        for(int e=rp[i];e<rp[i+1];e++) if(w[e]>mx) mx=w[e];
        int best=-1; double bw=-1.;
        for(int e=rp[i];e<rp[i+1];e++){
            int j=cj[e];
            if(j==i || agg[j]>=0 || !(w[e]>0.) || w[e]<GMG_BETA*mx) continue;
            if(w[e]>bw){ bw=w[e]; best=j; }
        }
        agg[i]=na; if(best>=0) agg[best]=na;
        na++;
    }
    return na;
}
/* contract a weighted graph by agg (na groups): CSR without self entries */
static int g_contract(int N, const int *rp, const int *cj, const double *w, const int *agg, int na,
                      int **orp, int **ocj, double **ow){
    int *cnt=calloc((size_t)na+1,sizeof(int)), *mem=malloc((size_t)(N>0?N:1)*sizeof(int));
    int *mark=malloc((size_t)(na>0?na:1)*sizeof(int)), *slot=malloc((size_t)(na>0?na:1)*sizeof(int));
    int cap=rp[N]>16?rp[N]:16, nnz=0;
    int *ncj=malloc((size_t)cap*sizeof(int)); double *nw=malloc((size_t)cap*sizeof(double));
    int *nrp=calloc((size_t)na+1,sizeof(int));
    if(!cnt||!mem||!mark||!slot||!ncj||!nw||!nrp){ free(cnt);free(mem);free(mark);free(slot);free(ncj);free(nw);free(nrp); return 0; }
    for(int i=0;i<N;i++) cnt[agg[i]+1]++;
    for(int a=0;a<na;a++) cnt[a+1]+=cnt[a];
    { int *pos=malloc((size_t)(na>0?na:1)*sizeof(int));
      if(!pos){ free(cnt);free(mem);free(mark);free(slot);free(ncj);free(nw);free(nrp); return 0; }
      for(int a=0;a<na;a++) pos[a]=cnt[a];
      for(int i=0;i<N;i++) mem[pos[agg[i]]++]=i;
      free(pos); }
    for(int a=0;a<na;a++) mark[a]=-1;
    for(int I=0;I<na;I++){
        for(int q=cnt[I];q<cnt[I+1];q++){
            int i=mem[q];
            for(int e=rp[i];e<rp[i+1];e++){
                int J=agg[cj[e]]; if(J==I) continue;
                if(mark[J]!=I){
                    if(nnz>=cap){ cap*=2;
                        int *t1=realloc(ncj,(size_t)cap*sizeof(int)); double *t2=realloc(nw,(size_t)cap*sizeof(double));
                        if(!t1||!t2){ free(t1?t1:ncj); free(t2?t2:nw); free(cnt);free(mem);free(mark);free(slot);free(nrp); return 0; }
                        ncj=t1; nw=t2; }
                    mark[J]=I; slot[J]=nnz; ncj[nnz]=J; nw[nnz]=0.; nnz++;
                }
                nw[slot[J]]+=w[e];
            }
        }
        nrp[I+1]=nnz;
    }
    free(cnt); free(mem); free(mark); free(slot);
    *orp=nrp; *ocj=ncj; *ow=nw;
    return 1;
}
static int g_aggregate(const Lev *f, int *agg){
    int N=(int)f->N;
    int *cur=malloc((size_t)N*sizeof(int)), *pa=malloc((size_t)N*sizeof(int));
    if(!cur||!pa){ free(cur); free(pa); return -1; }
    for(int i=0;i<N;i++) agg[i]=i;
    const int *rp=f->rp, *cj=f->cj; const double *w=f->wg;
    int *orp=NULL,*ocj=NULL; double *ow=NULL; int n=N, na=N;
    for(int pass=0;pass<GMG_PASSES && n>1;pass++){
        na=g_pair(n,rp,cj,w,pa);
        for(int i=0;i<N;i++) agg[i]=pa[agg[i]];
        if(pass==GMG_PASSES-1 || na==n) break;
        int *nrp,*ncj; double *nw;
        if(!g_contract(n,rp,cj,w,pa,na,&nrp,&ncj,&nw)){ na=-1; break; }
        free(orp); free(ocj); free(ow);
        orp=nrp; ocj=ncj; ow=nw; rp=orp; cj=ocj; w=ow; n=na;
    }
    free(orp); free(ocj); free(ow); free(cur); free(pa);
    return na;
}

/* coarse pattern of level c from the aggregation of f; member lists and
   the fine-entry -> coarse-entry map */
static int g_coarse_pattern(Lev *f, Lev *c){
    int N=(int)f->N, na=(int)c->N;
    f->mrp=calloc((size_t)na+1,sizeof(int)); f->mem=malloc((size_t)N*sizeof(int));
    f->emap=malloc((size_t)(f->nnz>0?f->nnz:1)*sizeof(int));
    int *mark=malloc((size_t)na*sizeof(int)), *slot=malloc((size_t)na*sizeof(int));
    if(!f->mrp||!f->mem||!f->emap||!mark||!slot){ free(mark); free(slot); return 0; }
    for(int i=0;i<N;i++) f->mrp[f->agg[i]+1]++;
    for(int a=0;a<na;a++) f->mrp[a+1]+=f->mrp[a];
    { int *pos=malloc((size_t)na*sizeof(int)); if(!pos){ free(mark); free(slot); return 0; }
      for(int a=0;a<na;a++) pos[a]=f->mrp[a];
      for(int i=0;i<N;i++) f->mem[pos[f->agg[i]]++]=i;
      free(pos); }
    int cap=f->nnz>16?f->nnz:16, nnz=0;
    int *cj=malloc((size_t)cap*sizeof(int));
    c->rp=calloc((size_t)na+1,sizeof(int));
    if(!cj||!c->rp){ free(mark); free(slot); free(cj); return 0; }
    for(int a=0;a<na;a++) mark[a]=-1;
    for(int I=0;I<na;I++){
        for(int q=f->mrp[I];q<f->mrp[I+1];q++){
            int i=f->mem[q];
            for(int e=f->rp[i];e<f->rp[i+1];e++){
                int J=f->agg[f->cj[e]];
                if(J==I){ f->emap[e]=-1; continue; }
                if(mark[J]!=I){
                    if(nnz>=cap){ cap*=2; int *t=realloc(cj,(size_t)cap*sizeof(int)); if(!t){ free(mark); free(slot); free(cj); return 0; } cj=t; }
                    mark[J]=I; slot[J]=nnz; cj[nnz++]=J;
                }
                f->emap[e]=slot[J];
            }
        }
        c->rp[I+1]=nnz;
    }
    free(mark); free(slot);
    c->cj=cj; c->nnz=nnz;
    c->wg=calloc((size_t)(nnz>0?nnz:1),sizeof(double));
    if(!c->wg) return 0;
    for(int i=0;i<N;i++) for(int e=f->rp[i];e<f->rp[i+1];e++) if(f->emap[e]>=0) c->wg[f->emap[e]]+=f->wg[e];
    c->rv=malloc((size_t)(nnz>0?nnz:1)*sizeof(int));
    if(!c->rv) return 0;
    for(int I=0;I<na;I++) for(int s=c->rp[I];s<c->rp[I+1];s++){
        int J=c->cj[s], r=-1;
        for(int t=c->rp[J];t<c->rp[J+1];t++) if(c->cj[t]==I){ r=t; break; }
        c->rv[s]=r;                              /* the pattern is symmetric: r >= 0 */
    }
    return 1;
}

static int gmg_alloc(Dev *d){
    mg_free(d);
    { const char *e=getenv("SEMISIM_GPASS"); if(e) GMG_PASSES=atoi(e);
      const char *b=getenv("SEMISIM_GBETA"); if(b) GMG_BETA=atof(b); }
    Lev *f=&d->lev[0];
    f->Nx=(int)NN(d); f->Ny=f->Nz=1; f->N=NN(d);
    f->rp=d->grp; f->cj=d->gcol; f->rv=d->grev; f->nnz=d->nE; f->va=d->gval; f->dg=d->A.c;
    f->ilu=d->ilu; f->r=calloc(f->N,8); f->pw=calloc(f->N,8);
    f->wg=malloc((size_t)(d->nE>0?d->nE:1)*sizeof(double));
    if(!f->r||!f->pw||!f->wg) return 0;
    for(int e=0;e<d->nE;e++) f->wg[e]= d->gh[e]>0. ? d->ga[e]/d->gh[e] : 0.;
    int nl=1;
    while(nl<MAXLEV){
        f=&d->lev[nl-1];
        if(f->N<=GMG_COARSE) break;
        f->agg=malloc(f->N*sizeof(int)); if(!f->agg) return 0;
        int na=g_aggregate(f,f->agg);
        if(na<0) return 0;
        if(na<1 || (double)na>0.7*(double)f->N){ free(f->agg); f->agg=NULL; break; }
        Lev *c=&d->lev[nl]; c->N=(size_t)na; c->Nx=na; c->Ny=c->Nz=1; c->own=1;
        if(!g_coarse_pattern(f,c)) return 0;
        c->va=calloc((size_t)(c->nnz>0?c->nnz:1),8); c->dg=calloc(c->N,8);
        c->ilu=calloc(c->N,8); c->x=calloc(c->N,8); c->b=calloc(c->N,8); c->r=calloc(c->N,8);
        c->pw=calloc(c->N,8); c->ph=calloc(c->N,8);
        if(!c->va||!c->dg||!c->ilu||!c->x||!c->b||!c->r||!c->pw||!c->ph) return 0;
        nl++;
    }
    Lev *z=&d->lev[nl-1];
    if(z->N<=4*GMG_COARSE){
        z->LU=calloc(z->N*z->N,8); z->piv=calloc(z->N,sizeof(int));
        if(!z->LU||!z->piv) return 0;
    }
    d->nlev=nl;
    if(d->verbose) fprintf(stderr,"semisim: graph multigrid %d levels, coarsest %zu unknowns\n",nl,z->N);
    return 1;
}

static void gilu_factor(const Dev *d, Lev *v){
    const size_t N=v->N; double *D=v->ilu;
    int nch=g_chunks(d,N);
    #pragma omp parallel for schedule(static,1) if(nch>1)
    for(int c=0;c<nch;c++){
        int k0=(int)(N*(size_t)c/nch), k1=(int)(N*(size_t)(c+1)/nch);
        for(int k=k0;k<k1;k++){
            double dk=v->dg[k];
            for(int e=v->rp[k];e<v->rp[k+1];e++){
                int j=v->cj[e];
                if(j<k && j>=k0 && v->va[e]!=0.) dk-=v->va[e]*v->va[v->rv[e]]/D[j];
            }
            if(!(dk>1e-8*fabs(v->dg[k]))) dk=v->dg[k];
            if(!(dk!=0.)) dk=1.;
            D[k]=dk;
        }
    }
}
static void gilu_apply(const Dev *d, const Lev *v, const double *r, double *z){
    const size_t N=v->N; const double *D=v->ilu;
    int nch=g_chunks(d,N);
    #pragma omp parallel for schedule(static,1) if(nch>1)
    for(int c=0;c<nch;c++){
        int k0=(int)(N*(size_t)c/nch), k1=(int)(N*(size_t)(c+1)/nch);
        for(int k=k0;k<k1;k++){
            double s=r[k];
            for(int e=v->rp[k];e<v->rp[k+1];e++){ int j=v->cj[e]; if(j<k && j>=k0) s-=v->va[e]*z[j]; }
            z[k]=s/D[k];
        }
        for(int k=k1-1;k>=k0;k--){
            double s=0.;
            for(int e=v->rp[k];e<v->rp[k+1];e++){ int j=v->cj[e]; if(j>k && j<k1) s+=v->va[e]*z[j]; }
            z[k]-=s/D[k];
        }
    }
}
static void gmatvec_L(const Dev *d, const Lev *v, const double *x, const double *b, double *y){
    const size_t N=v->N;
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        double s=v->dg[k]*x[k];
        for(int e=v->rp[k];e<v->rp[k+1];e++) s+=v->va[e]*x[v->cj[e]];
        y[k]= b ? b[k]-s : s;
    }
    (void)d;
}
/* rounding floor of the continuity system (see bicgstab_core) */
static double gfloor(const Dev *d, const double *b, const double *x){
    const Lev *v=&d->lev[0]; const size_t N=v->N; double fl=0.;
    #pragma omp parallel for reduction(+:fl) schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        double s=fabs(v->dg[k]*x[k])+fabs(b[k]);
        for(int e=v->rp[k];e<v->rp[k+1];e++) s+=fabs(v->va[e]*x[v->cj[e]]);
        fl+=s*s;
    }
    return fl;
}

static void g_galerkin(const Lev *f, Lev *c, const double *w){
    const int na=(int)c->N;
    #pragma omp parallel for schedule(static) if(c->N>PAR_MIN/8)
    for(int I=0;I<na;I++){
        for(int s=c->rp[I];s<c->rp[I+1];s++) c->va[s]=0.;
        double cc=0.;
        for(int q=f->mrp[I];q<f->mrp[I+1];q++){
            int i=f->mem[q]; double wi= w ? w[i] : 1.;
            cc+=wi*f->dg[i]*f->pw[i];
            for(int e=f->rp[i];e<f->rp[i+1];e++){
                double a=wi*f->va[e]*f->pw[f->cj[e]];
                if(f->emap[e]<0) cc+=a; else c->va[f->emap[e]]+=a;
            }
        }
        int zero = cc==0.;
        for(int s=c->rp[I];s<c->rp[I+1] && zero;s++) if(c->va[s]!=0.) zero=0;
        c->dg[I]= zero ? 1. : cc;               /* identity rows only: keep an identity row */
    }
}
static void g_restrict(const Lev *f, const Lev *c, const double *r, double *bc, const double *w){
    const int na=(int)c->N;
    #pragma omp parallel for schedule(static) if(c->N>PAR_MIN/8)
    for(int I=0;I<na;I++){
        double s=0.;
        for(int q=f->mrp[I];q<f->mrp[I+1];q++){ int i=f->mem[q]; s+= w ? w[i]*r[i] : r[i]; }
        bc[I]=s;
    }
}
static void g_prolong(const Lev *f, const double *xc, double *x){
    const size_t N=f->N;
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t i=0;i<N;i++) x[i]+=MG_ALPHA*f->pw[i]*xc[f->agg[i]];
}
static void g_dense_factor(Lev *v){
    int n=(int)v->N; double *M=v->LU;
    memset(M,0,(size_t)n*n*8);
    for(int k=0;k<n;k++){
        M[(size_t)k*n+k]=v->dg[k];
        for(int e=v->rp[k];e<v->rp[k+1];e++) M[(size_t)k*n+v->cj[e]]+=v->va[e];
    }
    for(int c=0;c<n;c++){
        int p=c; double mx=fabs(M[(size_t)c*n+c]);
        for(int r=c+1;r<n;r++) if(fabs(M[(size_t)r*n+c])>mx){ mx=fabs(M[(size_t)r*n+c]); p=r; }
        v->piv[c]=p;
        if(p!=c) for(int q=0;q<n;q++){ double t=M[(size_t)c*n+q]; M[(size_t)c*n+q]=M[(size_t)p*n+q]; M[(size_t)p*n+q]=t; }
        double dg=M[(size_t)c*n+c]; if(!(fabs(dg)>1e-300)) dg=M[(size_t)c*n+c]=1e-300;
        for(int r=c+1;r<n;r++){
            double f=M[(size_t)r*n+c]/dg; M[(size_t)r*n+c]=f;
            if(f!=0.) for(int q=c+1;q<n;q++) M[(size_t)r*n+q]-=f*M[(size_t)c*n+q];
        }
    }
}
static void g_weights(Dev *d){
    for(int L=0;L<d->nlev-1;L++){
        Lev *f=&d->lev[L], *c=&d->lev[L+1];
        if(d->mg_qsign==0){ for(size_t k=0;k<f->N;k++) f->pw[k]=1.; continue; }
        const int na=(int)c->N;
        #pragma omp parallel for schedule(static) if(c->N>PAR_MIN/8)
        for(int I=0;I<na;I++){
            double mx=PH_NONE;
            for(int q=f->mrp[I];q<f->mrp[I+1];q++){
                int i=f->mem[q]; double p= L==0 ? ph0(d,(size_t)i) : f->ph[i];
                if(p>mx) mx=p;
            }
            c->ph[I]=mx;
            for(int q=f->mrp[I];q<f->mrp[I+1];q++){
                int i=f->mem[q]; double p= L==0 ? ph0(d,(size_t)i) : f->ph[i];
                double e= (p<=0.5*PH_NONE) ? -800. : p-mx;
                f->pw[i]= e<-700. ? 0. : exp(e);
            }
        }
    }
}
static void gmg_setup(Dev *d){
    double t0=wtime();
    g_weights(d);
    for(int L=0;L<d->nlev-1;L++){
        g_galerkin(&d->lev[L],&d->lev[L+1], L==0 ? d->mg_w : NULL);
        gilu_factor(d,&d->lev[L]);
    }
    Lev *z=&d->lev[d->nlev-1];
    if(z->LU) g_dense_factor(z); else gilu_factor(d,z);
    d->t_mg+=wtime()-t0;
}
static void gvcycle(Dev *d, int L, const double *b, double *x){
    Lev *v=&d->lev[L];
    if(mg_off && L==0 && d->nlev>1){ gilu_apply(d,v,b,x); return; }
    if(L==d->nlev-1){
        if(v->LU) dense_solve(v,b,x);
        else {                                     /* coarsest too big for LU: ILU sweeps */
            gilu_apply(d,v,b,x);
            for(int s=0;s<4;s++){ gmatvec_L(d,v,x,b,v->r); gilu_apply(d,v,v->r,v->r);
                                  for(size_t k=0;k<v->N;k++) x[k]+=v->r[k]; }
        }
        return;
    }
    Lev *c=&d->lev[L+1]; size_t N=v->N;
    gilu_apply(d,v,b,x);
    gmatvec_L(d,v,x,b,v->r);
    g_restrict(v,c,v->r,c->b, L==0 ? d->mg_w : NULL);
    gvcycle(d,L+1,c->b,c->x);
    g_prolong(v,c->x,x);
    gmatvec_L(d,v,x,b,v->r);
    gilu_apply(d,v,v->r,v->r);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++) x[k]+=v->r[k];
}

/* ------------------------------------------------ graph geometry helpers
   Least-squares gradient of a node field over the node's semiconductor
   neighbours: g = M^-1 sum_e w_e (f_b - f_k) dr_e, M = sum_e w_e dr_e dr_e^T,
   w_e = 1/|dr|^2.  glsq holds M^-1 (symmetric: xx xy xz yy yz zz). */
static void g_lsq_setup(Dev *d){
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        double m[6]={0,0,0,0,0,0}; double *o=d->glsq+6*k;
        if(d->kind[k]==K_SEMI) for(int e=d->grp[k];e<d->grp[k+1];e++){
            size_t b=d->gcol[e]; if(d->kind[b]!=K_SEMI) continue;
            double dx=d->gpos[3*b]-d->gpos[3*k], dy=d->gpos[3*b+1]-d->gpos[3*k+1], dz=d->gpos[3*b+2]-d->gpos[3*k+2];
            double w=1./(dx*dx+dy*dy+dz*dz+1e-40);
            m[0]+=w*dx*dx; m[1]+=w*dx*dy; m[2]+=w*dx*dz; m[3]+=w*dy*dy; m[4]+=w*dy*dz; m[5]+=w*dz*dz;
        }
        double tr=m[0]+m[3]+m[5], eps=1e-6*tr+1e-30;
        m[0]+=eps; m[3]+=eps; m[5]+=eps;
        double a=m[0],b=m[1],c=m[2],dd=m[3],e=m[4],f=m[5];
        double A=dd*f-e*e, B=-(b*f-c*e), Cc=b*e-c*dd, D=a*f-c*c, E=-(a*e-b*c), F=a*dd-b*b;
        double det=a*A+b*B+c*Cc;
        if(!(fabs(det)>0.)){ for(int q=0;q<6;q++) o[q]=0.; continue; }
        o[0]=A/det; o[1]=B/det; o[2]=Cc/det; o[3]=D/det; o[4]=E/det; o[5]=F/det;
    }
}
static inline void g_grad(const Dev *d, const double *fld, size_t k, double *g){
    double r[3]={0,0,0};
    if(d->kind[k]==K_SEMI) for(int e=d->grp[k];e<d->grp[k+1];e++){
        size_t b=d->gcol[e]; if(d->kind[b]!=K_SEMI) continue;
        double dx=d->gpos[3*b]-d->gpos[3*k], dy=d->gpos[3*b+1]-d->gpos[3*k+1], dz=d->gpos[3*b+2]-d->gpos[3*k+2];
        double w=(fld[b]-fld[k])/(dx*dx+dy*dy+dz*dz+1e-40);
        r[0]+=w*dx; r[1]+=w*dy; r[2]+=w*dz;
    }
    const double *m=d->glsq+6*k;
    g[0]=m[0]*r[0]+m[1]*r[1]+m[2]*r[2];
    g[1]=m[1]*r[0]+m[3]*r[1]+m[4]*r[2];
    g[2]=m[2]*r[0]+m[4]*r[1]+m[5]*r[2];
}

/* ============================================================ Poisson
   F_k = sum_f eps_f A_f (phi_nb-phi_k)/h + V_k q (p-n+Nd-Na+Nf)
         + A_g [C_ox (V_G' - phi_k) + q_ox]            (oxide gates)
   Newton: (-dF/dphi) dphi = F, SPD; PCG + IC(0); backtracking on ||F/diag||. */
static inline int is_dir(const Dev *d, size_t k){ return d->con_mask[k]!=0; }
/* unstructured-mesh variants (part3g.c) */
static void poisson_core_g(Dev *d, const double *phi, int build, double *F);
static void upd_mob_ref_g(Dev *d, const double *qn_, const double *qp_);
static void mob_display_g(Dev *d);
static void cont_assemble_g(Dev *d, int hole);
static int  float_fix_g(Dev *d, int hole);
static void terminal_flux_g(Dev *d, double *I, double *S);
static void upd_J_g(Dev *d);
static void region_shift_g(Dev *d, const double *dV);

static double poisson_eval(Dev *d, const double *phi, const double *Vapp,
                           int build, double *F){
    const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz; Sten *A=&d->A; size_t N=NN(d);
    if(d->unstr) poisson_core_g(d,phi,build,F);
    else
    #pragma omp parallel for collapse(2) schedule(static) if(N>PAR_MIN)
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        if(is_dir(d,k)){
            F[k]=0.;
            if(build){ A->c[k]=1.; A->xm[k]=A->xp[k]=A->ym[k]=A->yp[k]=A->zm[k]=A->zp[k]=0.;
                       d->sc[k]=1.; }
            continue;
        }
        double ek=d->eps[k], pk=phi[k], f=0., dg=0., off[6];
        for(int dir=0;dir<6;dir++){
            size_t kb; double h,ar; off[dir]=0.;
            if(!nbr(d,i,j,l,dir,&kb,&h,&ar)) continue;
            double eb=d->eps[kb], cf=(2.*ek*eb/(ek+eb))*ar/h;
            f+=cf*(phi[kb]-pk); dg+=cf;
            if(!is_dir(d,kb)) off[dir]=-cf;
        }
        double vol=volume(d,i,j,l);
        if(d->kind[k]==K_SEMI){
            const Mat *m=&d->mats[d->mat_idx[k]];
            double n=nK(d,k,pk,d->phi_n[k]), p=pK(d,k,pk,d->phi_p[k]);
            f +=vol*Q*(p-n+d->Nd[k]-d->Na[k]+d->Nf[k]);
            dg+=vol*Q*(n+p)/m->VT;
        } else f+=vol*Q*d->Nf[k];
        F[k]=f;
        if(build){
            A->c[k]=dg; A->xm[k]=off[0]; A->xp[k]=off[1]; A->ym[k]=off[2];
            A->yp[k]=off[3]; A->zm[k]=off[4]; A->zp[k]=off[5]; d->sc[k]=dg;
        }
    }
    for(int r=0;r<d->n_rob;r++){
        size_t k=d->rob_k[r]; int c=d->rob_c[r]; const Con *cn=&d->cons[c];
        double Cox=EPS_OX/(cn->tox>0.?cn->tox:10e-9);
        double pm = cn->phi_m>0. ? cn->phi_m : d->Cref;
        double VG = Vapp[c]+d->Cref-pm;
        F[k]+=d->rob_A[r]*(Cox*(VG-phi[k])+cn->qox);
        if(build){ A->c[k]+=Cox*d->rob_A[r]; d->sc[k]=A->c[k]; }
    }
    double ss=0.;
    #pragma omp parallel for reduction(+:ss) schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){ double r=F[k]/d->sc[k]; ss+=r*r; }
    return sqrt(ss);
}

static int poisson_newton(Dev *d, const double *Vapp, double ntol, int maxit){
    size_t N=NN(d);
    double *dx=d->sol, *tr=d->wk[4], *Ft=d->wk[5];
    for(int it=0;it<maxit;it++){
        d->newton_its++;
        double f0=poisson_eval(d,d->phi,Vapp,1,d->rhs);
        if(!(f0>0.)) return f0==0.;
        memset(dx,0,N*sizeof(double));
        d->mg_w=NULL; d->mg_qsign=0;
        pcg(d,d->rhs,dx,d->max_iter_p>0?d->max_iter_p:2000,d->pcg_rtol);
        double dm=0.;
        #pragma omp parallel for reduction(max:dm) schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){ double a=fabs(dx[k]); if(a>dm) dm=a; }
        if(dm!=dm) return 0;
        if(dm<ntol){
            #pragma omp parallel for schedule(static) if(N>PAR_MIN)
            for(size_t k=0;k<N;k++) d->phi[k]+=dx[k];
            return 1;
        }
        double lam=1.;
        for(int ls=0;ls<40;ls++){
            #pragma omp parallel for schedule(static) if(N>PAR_MIN)
            for(size_t k=0;k<N;k++) tr[k]=d->phi[k]+lam*dx[k];
            double f1=poisson_eval(d,tr,Vapp,0,Ft);
            if(f1<=(1.-1e-4*lam)*f0) break;
            lam*=0.5;
        }
        memcpy(d->phi,tr,N*sizeof(double));
        if(d->verbose>4) fprintf(stderr,"      newton it=%d f0=%.3e dm=%.3e lam=%.3e\n",it,f0,dm,lam);
        if(lam*dm<ntol) return 1;
    }
    return 0;
}

/* ======================================================= carriers */
static void rebuild_np(Dev *d){
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        if(d->kind[k]!=K_SEMI){ d->n[k]=d->p[k]=0.; continue; }
        d->n[k]=nK(d,k,d->phi[k],d->phi_n[k]);
        d->p[k]=pK(d,k,d->phi[k],d->phi_p[k]);
    }
}

/* ================================================ fluxes and currents
   Current (A) through face a->b of area ar, split into electron/hole parts;
   S* return the magnitude scale (sum of the two SG terms). */
static inline int face_axis(const Dev *d, size_t ka, size_t kb){
    size_t df = ka<kb ? kb-ka : ka-kb;
    return df==1 ? 0 : (df==(size_t)d->Nx ? 1 : 2);
}
static inline void face_flux(const Dev *d, size_t ka, size_t kb, double h, double ar,
                             double *Jn, double *Jp, double *Sn, double *Sp){
    const Mat *ma=&d->mats[d->mat_idx[ka]], *mb=&d->mats[d->mat_idx[kb]];
    double VT=ma->VT, dps=(d->phi[kb]-d->phi[ka])/VT; (void)mb;
    double Dn=dps+dlnc(d,ka,kb), Dp=dps-dlnv(d,ka,kb);
    int ax=face_axis(d,ka,kb); size_t kl= ka<kb ? ka : kb;
    double gn=Q*ar*(double)d->fmu[ax][kl]*VT/h;
    double gp=Q*ar*(double)d->fmu[3+ax][kl]*VT/h;
    double a1=gn*d->n[kb]*B(Dn), a2=gn*d->n[ka]*B(-Dn);
    double b1=gp*d->p[ka]*B(Dp), b2=gp*d->p[kb]*B(-Dp);
    *Jn=a1-a2; *Jp=b1-b2; *Sn=a1+a2; *Sp=b1+b2;
}

/* Mobility. Node arrays mu_n/mu_p hold the doping-dependent low-field value
   during a solve. Each SG edge gets its own mobility, saturated
   (Caughey-Thomas) by a driving force derived from carrier heating: the
   Joule power per carrier is J.E, and J = -q mu c grad(phi_c), so the local
   equivalent field is F^2 = (-grad phi_c).E, i.e. along an edge
       F = sqrt(max(0, d phi_c * d psi)) / h        (c = n or p).
   Drift (d phi_c = d psi): F = |E|. Forward-biased or high-low junctions,
   where drift and diffusion balance (d phi_c ~ 0): F ~ 0, so built-in fields
   do not throttle injection. Pure diffusion (d psi ~ 0): F ~ 0, so diffusing
   minority carriers are not "heated". Diffusion against a retarding field
   (product < 0): F = 0. This is the local limit of the energy-balance model. */
/* Driving field of an edge for velocity saturation: carriers are heated by
   the field only as far as they actually drift along it, so
       F = min(|d phi_c|, |d psi|)/h   if d phi_c and d psi have equal sign,
       F = 0                           otherwise.
   Drift (d phi_c = d psi): F = |E|. Near-equilibrium junctions (d phi_c ~ 0,
   drift balanced by diffusion): F ~ 0, so built-in fields do not throttle
   injection. Pure diffusion (d psi ~ 0): F ~ 0 (no heating of diffusing
   minority carriers). Against a retarding field: F = 0. Piecewise linear in
   both potentials, so it stays well conditioned inside the Gummel loop. */
/* Velocity-saturation driving field on an edge: the harmonic mean of the
   quasi-Fermi and electrostatic potential drops when both push the carrier
   the same way, else 0.  Equals the field in drift-dominated channels
   (dq = dpsi), -> 0 at equilibrium (dq -> 0, no built-in-field artefact), and
   is smooth where dq ~ dpsi (a min() switches branch there between Gummel
   iterations and stalls the iteration in saturated channels). */
static inline double heat_field(double dq, double dps, double h){
    if(dq*dps<=0.) return 0.;
    double a=fabs(dq), b=fabs(dps);
    return 2.*a*b/((a+b)*h);
}
static void mob0_init(Dev *d){
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        if(d->kind[k]!=K_SEMI){ d->mu_n[k]=d->mu_p[k]=0.; continue; }
        const Mat *m=MK(d,k); double Nt=d->Nd[k]+d->Na[k];
        d->mu_n[k]=mob_low(m,Nt,0); d->mu_p[k]=mob_low(m,Nt,1);
    }
}
/* Face mobilities from phi (current iterate) and the quasi-Fermi REFERENCE
   potentials qn/qp. During a Gummel solve the reference is frozen, so the
   Gummel map is a pure function of phi (needed by Anderson acceleration and
   free of the neutrally-stable lagged-mobility mode of saturated channels);
   an outer loop (mob_polish) updates the reference to self-consistency. */
/* Lombardi surface mobility (Lombardi et al. 1988, Sentaurus parameters,
   cm units):  1/mu = 1/mu_b + D/mu_ac + D/mu_sr,  D = exp(-d/10 nm),
     mu_ac = B/F + C N^lambda / (F^(1/3) (T/300))   (surface phonons)
     mu_sr = delta / F^2                            (surface roughness)
   F = field normal to the edge (V/cm): the potential gradient along the two
   axes perpendicular to it, taken inside the semiconductor only (one-sided
   at an interface - the field in the oxide is not the channel field). */
#define LB_LCRIT 1e-8
static inline double grad_semi(const Dev *d, int i, int j, int l, int ax){
    size_t k=IDX(i,j,l,d->Nx,d->Ny), km=0,kp=0; double hm=1.,hp=1.,ar;
    int a = nbr(d,i,j,l,2*ax,&km,&hm,&ar) && d->kind[km]==K_SEMI;
    int b = nbr(d,i,j,l,2*ax+1,&kp,&hp,&ar) && d->kind[kp]==K_SEMI;
    if(a&&b) return (d->phi[kp]-d->phi[km])/(hm+hp);
    if(b) return (d->phi[kp]-d->phi[k])/hp;
    if(a) return (d->phi[k]-d->phi[km])/hm;
    return 0.;
}
static inline double fnorm2(const Dev *d, int i, int j, int l, int ax){
    double g1=grad_semi(d,i,j,l,(ax+1)%3), g2=grad_semi(d,i,j,l,(ax+2)%3);
    return g1*g1+g2*g2;
}
static inline double lombardi(const Mat *m, double mub, double Fv, double Nt, double Dw, int hole){
    double F=Fv*1e-2; if(F<1.) F=1.;                       /* V/cm */
    double N=Nt*1e-6; if(N<1.) N=1.;                       /* cm^-3 */
    double Bc=hole?m->lbB_p:m->lbB_n, Cc=hole?m->lbC_p:m->lbC_n;
    double Lc=hole?m->lbL_p:m->lbL_n, Dc=hole?m->lbD_p:m->lbD_n;
    double muac=Bc/F+Cc*pow(N,Lc)/(cbrt(F)*(m->T/300.));
    double musr=Dc/(F*F);
    return 1./(1./mub+Dw*1e4*(1./muac+1./musr));           /* cm^2/Vs -> m^2/Vs */
}
static void upd_mob_ref(Dev *d, const double *qn_, const double *qp_){
    const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz;
    const int sat = (d->models&4)!=0, lb = (d->models&16)!=0;
    #pragma omp parallel for collapse(2) schedule(static) if(NN(d)>PAR_MIN)
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        for(int ax=0;ax<3;ax++){
            size_t kp; double h,ar;
            if(d->kind[k]!=K_SEMI || !nbr(d,i,j,l,2*ax+1,&kp,&h,&ar) || d->kind[kp]!=K_SEMI){
                d->fmu[ax][k]=0.; d->fmu[3+ax][k]=0.; continue;
            }
            const Mat *m=MK(d,k);
            double mn=0.5*(d->mu_n[k]+d->mu_n[kp]), mp=0.5*(d->mu_p[k]+d->mu_p[kp]);
            if(lb && m->lb_on){
                double Dw=exp(-0.5*(d->qd[k]+d->qd[kp])/LB_LCRIT);
                if(Dw>1e-8){
                    int ip=i+(ax==0), jp=j+(ax==1), lp=l+(ax==2);
                    double F=sqrt(0.5*(fnorm2(d,i,j,l,ax)+fnorm2(d,ip,jp,lp,ax)));
                    double Nt=0.5*(d->Nd[k]+d->Na[k]+d->Nd[kp]+d->Na[kp]);
                    mn=lombardi(m,mn,F,Nt,Dw,0); mp=lombardi(m,mp,F,Nt,Dw,1);
                }
            }
            if(sat){
                double dps=d->phi[kp]-d->phi[k];
                double Fn=heat_field(qn_[kp]-qn_[k],dps,h), Fp=heat_field(qp_[kp]-qp_[k],dps,h);
                mn=mob_hf(mn,Fn,m->vsat_n,m->beta_n); mp=mob_hf(mp,Fp,m->vsat_p,m->beta_p);
            }
            d->fmu[ax][k]=(mu_t)mn; d->fmu[3+ax][k]=(mu_t)mp;
        }
    }
}
static void upd_mob(Dev *d){
    if(d->unstr){ upd_mob_ref_g(d,d->phi_n,d->phi_p); return; }
    if(d->verbose>2 && (d->models&4)){
        /* diagnostic: largest relative change of an edge mobility vs last call */
        const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz; size_t N=NN(d);
        mu_t *old=(mu_t*)malloc(6*N*sizeof(mu_t));
        for(int q=0;q<6;q++) memcpy(old+q*N,d->fmu[q],N*sizeof(mu_t));
        upd_mob_ref(d,d->phi_n,d->phi_p);
        double mx=0; size_t km=0; int qm=0;
        for(int q=0;q<6;q++) for(size_t k=0;k<N;k++){
            double a=d->fmu[q][k], b=old[q*N+k]; if(a<=0) continue;
            double w = q<3 ? d->n[k] : d->p[k];
            if(w<1e-6*(d->n[k]+d->p[k]+fabs(d->Nd[k]-d->Na[k]))) continue;
            double r=fabs(a-b)/a; if(r>mx){ mx=r; km=k; qm=q; }
        }
        int i,j,l; unidx(d,km,&i,&j,&l);
        size_t kp = km + (qm%3==0 ? 1 : (qm%3==1 ? (size_t)Nx : (size_t)Nx*Ny));
        fprintf(stderr,"    mob change %.3e at (%d,%d,%d) %s-edge %c: dphi=%.4e dqf=%.4e n=%.2e mu=%.1f\n",
            mx,i,j,l,(qm%3==0?"x":(qm%3==1?"y":"z")),qm<3?'n':'p',
            kp<N? d->phi[kp]-d->phi[km]:0., kp<N? (qm<3? d->phi_n[kp]-d->phi_n[km] : d->phi_p[kp]-d->phi_p[km]):0.,
            d->n[km], d->fmu[qm][km]*1e4);
        free(old); (void)Nz;
        return;
    }
    upd_mob_ref(d,d->phi_n,d->phi_p);
}
/* display mobility: mean of the adjacent edge mobilities (call after solve) */
static void mob_display(Dev *d){
    if(d->unstr){ mob_display_g(d); return; }
    const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz;
    #pragma omp parallel for collapse(2) schedule(static) if(NN(d)>PAR_MIN)
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        if(d->kind[k]!=K_SEMI) continue;
        double sn=0.,sp=0.; int c=0;
        for(int dir=0;dir<6;dir++){
            size_t kb; double h,ar;
            if(!nbr(d,i,j,l,dir,&kb,&h,&ar) || d->kind[kb]!=K_SEMI) continue;
            size_t kl= k<kb?k:kb; int ax=dir/2;
            sn+=d->fmu[ax][kl]; sp+=d->fmu[3+ax][kl]; c++;
        }
        if(c){ d->wk[0][k]=sn/c; d->wk[1][k]=sp/c; }
        else { d->wk[0][k]=d->mu_n[k]; d->wk[1][k]=d->mu_p[k]; }
    }
    size_t N=NN(d);
    for(size_t k=0;k<N;k++) if(d->kind[k]==K_SEMI){ d->mu_n[k]=d->wk[0][k]; d->mu_p[k]=d->wk[1][k]; }
}

/* Scaled SG continuity. Electrons (hole=0), holes (hole=1).
   Row a (unscaled):  sum_b [self_ab c_a - nb_ab c_b] + V_a dR c_a = V_a (dR c_a0 - R_a)
     electrons: D = dphi/VT + d(lnc), self = g B(-D), nb = g B(D)
     holes:     D = dphi/VT - d(lnv), self = g B(D),  nb = g B(-D)
   Unknown u = c/s with s = current iterate; row scaled by 1/(diag s_a). */
static void cont_assemble(Dev *d, int hole){
    if(d->unstr){ cont_assemble_g(d,hole); return; }
    const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz; Sten *A=&d->A;
    const double *car = hole? d->p : d->n;
    #pragma omp parallel for collapse(2) schedule(static) if(NN(d)>PAR_MIN)
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        if(d->kind[k]!=K_SEMI || d->con_mask[k]){
            A->c[k]=1.; A->xm[k]=A->xp[k]=A->ym[k]=A->yp[k]=A->zm[k]=A->zp[k]=0.;
            d->rhs[k]=(d->kind[k]==K_SEMI)?1.:0.; d->sc[k]=(d->kind[k]==K_SEMI)?car[k]:0.;
            d->dsc[k]=0.;   /* identity rows carry no flux: excluded from coarse levels */
            continue;
        }
        const Mat *ma=&d->mats[d->mat_idx[k]];
        double VT=ma->VT, sk=fmax(car[k],NFLOOR), dg=0., b=0., off[6];
        for(int dir=0;dir<6;dir++){
            size_t kb; double h,ar; off[dir]=0.;
            if(!nbr(d,i,j,l,dir,&kb,&h,&ar) || d->kind[kb]!=K_SEMI) continue;
            double D=(d->phi[kb]-d->phi[k])/VT + (hole ? -dlnv(d,k,kb) : dlnc(d,k,kb));
            size_t kl= k<kb ? k : kb;
            double g=ar*(double)d->fmu[(hole?3:0)+dir/2][kl]*VT/h;
            double self = hole ? g*B(D)  : g*B(-D);
            double nb   = hole ? g*B(-D) : g*B(D);
            dg+=self;
            if(d->con_mask[kb]) b+=nb*car[kb];
            else off[dir]=-nb*fmax(car[kb],NFLOOR);
        }
        double R,dRn,dRp; recombK(d,k,&R,&dRn,&dRp);
        double dR = hole?dRp:dRn; if(dR<0.) dR=0.;
        double vol=volume(d,i,j,l);
        dg+=vol*dR; b+=vol*(dR*car[k]-R);
        dg+=d->tk[k]; b+=d->tg[k];                 /* gate tunnelling: sink tk*c - tg */
        double w=1./(dg*sk);
        A->c[k]=1.; A->xm[k]=off[0]*w; A->xp[k]=off[1]*w; A->ym[k]=off[2]*w;
        A->yp[k]=off[3]*w; A->zm[k]=off[4]*w; A->zp[k]=off[5]*w;
        d->rhs[k]=b*w; d->sc[k]=sk; d->dsc[k]=dg*sk;
    }
    d->mg_w=d->dsc; d->mg_qsign= hole ? -1 : 1;
}

/* c_new = s*u -> quasi-Fermi potential (clamped to the terminal range).
   Stores the change in dq; returns the max weighted |change| (weight 1 for a
   carrier >= 1e-3 of n+p+|Nd-Na|, linear down to 1e-3 at 1e-6, 0 below).
   Trace carriers are left to the separate terminal-current test. */
static double cont_update(Dev *d, int hole, double *dq, double lo, double hi){
    size_t N=NN(d); double dmax=0.;
    #pragma omp parallel for reduction(max:dmax) schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        if(d->kind[k]!=K_SEMI || d->con_mask[k]){ dq[k]=0.; continue; }
        const Mat *m=&d->mats[d->mat_idx[k]];
        double u=d->sol[k];
        if(!(u>1e-30)) u=1e-30; else if(u>1e30) u=1e30;
        double c=d->sc[k]*u, qn;
        if(hole){ qn=d->phi[k]+m->VT*(log(c)-lnvK(d,k)); }
        else    { qn=d->phi[k]-m->VT*(log(c)-lncK(d,k)); }
        if(qn<lo) qn=lo; else if(qn>hi) qn=hi;
        double *qf = hole? d->phi_p : d->phi_n;
        double ch=qn-qf[k]; qf[k]=qn; dq[k]=ch;
        double cn = hole ? pK(d,k,d->phi[k],qn) : nK(d,k,d->phi[k],qn);
        if(hole) d->p[k]=cn; else d->n[k]=cn;
        /* convergence weight: 1 where the carrier is >= 1e-3 of the local
           charge scale, falling linearly to 1e-3 at 1e-6 and 0 below: a
           trace carrier barely touches the electrostatics, while its
           current contribution is tested by the terminal-current check. */
        double scale=d->n[k]+d->p[k]+fabs(d->Nd[k]-d->Na[k]);
        double r = scale>0. ? cn/scale : 0.;
        if(r>1e-6){ double w = r>=1e-3 ? 1. : r*1e3, a=w*fabs(ch); if(a>dmax) dmax=a; }
    }
    if(d->verbose>2){
        size_t kx=0; double mx=-1;
        for(size_t k=0;k<N;k++){ if(d->kind[k]!=K_SEMI||d->con_mask[k]) continue;
            const Mat *m=&d->mats[d->mat_idx[k]]; double cc=hole?d->p[k]:d->n[k];
            if(cc>1e-6*(d->n[k]+d->p[k]+fabs(d->Nd[k]-d->Na[k])) && fabs(dq[k])>mx){ mx=fabs(dq[k]); kx=k; } (void)m; }
        int i,j,l; unidx(d,kx,&i,&j,&l);
        fprintf(stderr,"    %s max dq=%.3e at (%d,%d,%d) n=%.3e p=%.3e phi=%.4f phin=%.4f phip=%.4f\n",
            hole?"p":"n",mx,i,j,l,d->n[kx],d->p[kx],d->phi[kx],d->phi_n[kx],d->phi_p[kx]);
    }
    return dmax;
}

/* Floating majority-carrier regions.  A connected region where carrier c
   is the majority and that no carrier contact reaches (the n- base of a
   blocking thyristor or of an IGBT with its gate off, an inversion layer
   without source/drain) talks to the rest of the device only through
   exponentially small junction currents.  Its uniform quasi-Fermi shift is
   then the near-null mode of the continuity system and the iterative solve
   leaves it at round-off noise (the Gummel loop never settles).  After each
   solve every such block is rescaled uniformly, c -> lambda*c, so that its
   exact global balance holds - net outflow through its boundary plus
   recombination inside = 0 - a scalar equation whose terms are all
   evaluated directly (no cancellation), solved by Newton in ln(lambda). */
static int float_fix(Dev *d, int hole){
    if(d->unstr) return float_fix_g(d,hole);
    size_t N=NN(d);
    int *lab=(int*)d->wk[0]; size_t *q=(size_t*)d->wk[1];
    const double *oth = hole ? d->n : d->p;
    double *u=d->sol, *s=d->sc;
    /* majority AND quasi-neutral: c >= other, c >= |Nd-Na|/2 and well above
       n_i - depleted spots, where 'majority' flips on noise, never qualify */
    #define MAJ(k) (d->kind[k]==K_SEMI && !d->con_mask[k] && s[k]*u[k]>=oth[k] \
                    && s[k]*u[k]>=0.5*fabs(d->Nd[k]-d->Na[k]) && s[k]*u[k]>10.*niK(d,k))
    for(size_t k=0;k<N;k++) lab[k]=0;
    int nfix=0, blk=0;
    for(size_t k0=0;k0<N;k0++){
        if(lab[k0] || !MAJ(k0)) continue;
        blk++; size_t head=0,tail=0; q[tail++]=k0; lab[k0]=blk; int anch=0;
        while(head<tail){
            size_t k=q[head++]; int i,j,l; unidx(d,k,&i,&j,&l);
            for(int dir=0;dir<6;dir++){
                size_t kb; double h,ar;
                if(!nbr(d,i,j,l,dir,&kb,&h,&ar) || d->kind[kb]!=K_SEMI) continue;
                if(d->con_mask[kb]){ anch=1; continue; }
                if(!lab[kb] && MAJ(kb)){ lab[kb]=blk; q[tail++]=kb; }
            }
        }
        if(anch) continue;
        /* balance of the block: G(lam) = lam*Aout - Bin + sum vol*R(lam*c, other) */
        double Aout=0., Bin=0.;
        for(size_t t=0;t<tail;t++){
            size_t k=q[t]; int i,j,l; unidx(d,k,&i,&j,&l);
            const Mat *ma=&d->mats[d->mat_idx[k]]; double VT=ma->VT;
            for(int dir=0;dir<6;dir++){
                size_t kb; double h,ar;
                if(!nbr(d,i,j,l,dir,&kb,&h,&ar) || d->kind[kb]!=K_SEMI || lab[kb]==blk) continue;
                double D=(d->phi[kb]-d->phi[k])/VT + (hole ? -dlnv(d,k,kb) : dlnc(d,k,kb));
                size_t kl= k<kb ? k : kb;
                double g=ar*(double)d->fmu[(hole?3:0)+dir/2][kl]*VT/h;
                double self = hole ? g*B(D) : g*B(-D), nb = hole ? g*B(-D) : g*B(D);
                double cb = d->con_mask[kb] ? (hole?d->p[kb]:d->n[kb]) : s[kb]*u[kb];
                Aout+=self*s[k]*u[k]; Bin+=nb*cb;
            }
        }
        double sl=0.;
        for(int it=0;it<60;it++){
            double lam=exp(sl), G=lam*Aout-Bin, dG=lam*Aout;
            for(size_t t=0;t<tail;t++){
                size_t k=q[t]; int i,j,l; unidx(d,k,&i,&j,&l);
                double c=lam*s[k]*u[k], R,dRn,dRp;
                recomb(MK(d,k),niK(d,k),tauf(d,k), hole?oth[k]:c, hole?c:oth[k], &R,&dRn,&dRp);
                double vol=volume(d,i,j,l);
                G+=vol*R; dG+=vol*(hole?dRp:dRn)*c;
                G+=d->tk[k]*c-d->tg[k]; dG+=d->tk[k]*c;
            }
            if(!(dG>0.)) break;
            double ds=-G/dG; if(ds>2.) ds=2.; else if(ds<-2.) ds=-2.;
            sl+=ds; if(fabs(ds)<1e-13) break;
        }
        if(sl==sl && fabs(sl)>1e-15 && fabs(sl)<200.){
            double lam=exp(sl);
            for(size_t t=0;t<tail;t++) u[q[t]]*=lam;
            nfix++;
            if(d->verbose>2) fprintf(stderr,"    float_fix(%c) block of %zu nodes: shift %.3e V\n",hole?'p':'n',tail,
                                     (hole?-1.:1.)*d->mats[d->mat_idx[q[0]]].VT*sl);
        }
    }
    #undef MAJ
    return nfix;
}

/* ================================================== gate tunnelling
   Tsu-Esaki current with WKB transmission through the layer stack of each
   path (semiconductor node -> gate metal):
     F = 4 pi m_s kT/h^3  Int T(E) ln(1+exp((E_F-E)/kT)) dE     (1/m^2 s)
   evaluated separately for the supply of the semiconductor (s->g, its
   quasi-Fermi level) and of the metal (g->s): electrons of the conduction
   band (E >= Ec'), holes of the valence band (E <= Ev', mirrored).  The
   field follows from the series capacitance between the node and the
   metal, D = (psi_m - psi_k)/sum(t_i/eps_i), psi linear in every layer;
   Ec_ox = -psi - chi_ox, Ev_ox = Ec_ox - Eg_ox (vacuum level = -psi).
   Direct and Fowler-Nordheim tunnelling are the same integral (the barrier
   turns triangular once the oxide drop exceeds its height); over-barrier
   emission is included (T = 1).
   Supply.  Carriers pressed against the interface by the field (inversion
   or accumulation) supply the current through the rate at which they hit
   it, i.e. through their pressure on the wall.  For Boltzmann carriers in
   quasi-equilibrium along the normal the contact theorem gives it exactly,
       kT n(0) = kT n(z1) + q Int_0^z1 n F dz ,
   summed over the column of nodes inward from the interface.  This fixes
   the effective band edge of the supply, n_eff = Nc exp((E_F-Ec')/kT),
   independently of the mesh and of the quantum correction (MLDA empties
   the first nm but leaves the pressure - sheet charge x field - intact).
   A carrier pushed away from the interface uses its local density,
   extrapolated to the interface.  Both directions use the same Ec', so the
   net current vanishes exactly in equilibrium.
   Continuity: the escape is distributed over the column in proportion to
   each node's share of the pressure (exactly linear in the densities for
   Boltzmann supply), the injection G enters at the interface node. */
static inline double softplus(double x){ return x>35. ? x+log1p(exp(-x)) : log1p(exp(x)); }
/* 2 Int kappa dx over a layer of length L whose barrier above the carrier
   energy runs linearly from u0 to u1 (eV); kappa = sqrt(2 m q u)/hbar */
static inline double wkb_seg(double m, double L, double u0, double u1){
    double a=u0>0.?u0:0., b=u1>0.?u1:0.;
    if(a==0. && b==0.) return 0.;
    double c=2.*sqrt(2.*m*M0*Q)/HBAR, du=u1-u0, I;
    if(fabs(du)<1e-9*(a+b)) I=L*sqrt(0.5*(a+b));
    else I=L*(2./3.)*(b*sqrt(b)-a*sqrt(a))/du;
    return c*I;
}
#define TUN_NE   64       /* Simpson intervals per energy integral */
#define TUN_EMAX 3.0      /* integration range beyond the band edge (eV) */
typedef struct { int nl; const Mat *m[TUN_MAXL]; double t[TUN_MAXL], p0[TUN_MAXL], dp[TUN_MAXL];
                 double D, psis; } TPath;
static double tun_T(const TPath *P, double E, int hole){
    double s=0.;
    for(int i=0;i<P->nl;i++){
        const Mat *m=P->m[i]; if(!m->insulator) continue;
        double mt = hole? m->mt_h : m->mt_e; if(!(mt>0.)) mt=0.5;
        double ea=-P->p0[i]-m->chi, eb=-(P->p0[i]+P->dp[i])-m->chi;   /* Ec_ox at both ends */
        double u0,u1;
        if(hole){ u0=E-(ea-m->Eg); u1=E-(eb-m->Eg); } else { u0=ea-E; u1=eb-E; }
        s+=wkb_seg(mt,P->t[i],u0,u1);
        if(s>745.) return 0.;
    }
    return exp(-s);
}
/* Int T(E) ln(1+exp(+-(EF-E)/kT)) dE (eV): electrons upwards from E0 = Ec',
   holes downwards from E0 = Ev' */
static double tun_int(const TPath *P, double E0, double EF, double kT, int hole){
    double sg = hole ? -1. : 1.;
    double top = sg*(EF-E0); if(top<0.) top=0.;
    double span = top+25.*kT; if(span>TUN_EMAX) span=TUN_EMAX;
    double h=span/TUN_NE, s=0.;
    for(int q=0;q<=TUN_NE;q++){
        double E=E0+sg*q*h;
        double f=tun_T(P,E,hole)*softplus(sg*(EF-E)/kT);
        s+=f*(q==0||q==TUN_NE ? 1. : (q&1 ? 4. : 2.));
    }
    return s*h/3.;
}
/* potentials along a stack (semiconductor half cell first, if any) between a
   node at psik and a metal at voltage Vg with work function pm */
static void tun_stack(const Dev *d, TPath *P, int nl, const int *mi, const double *t,
                      double psik, double Vg, double pm){
    P->nl=nl; double R=0.;
    for(int i=0;i<nl;i++){ P->m[i]=&d->mats[mi[i]]; P->t[i]=t[i]; R+=t[i]/(P->m[i]->eps_r*EPS0); }
    if(!(pm>0.)) pm=d->Cref;
    double psim=Vg+d->Cref-pm, ps=psik;
    P->D = R>0. ? (psim-psik)/R : 0.;
    P->psis=psik;
    for(int i=0;i<nl;i++){
        P->p0[i]=ps; P->dp[i]=P->D*t[i]/(P->m[i]->eps_r*EPS0); ps+=P->dp[i];
        if(i==0 && !P->m[0]->insulator) P->psis=ps;            /* potential at the interface */
    }
}
/* fluxes (1/m^2 s) for effective band edges Ec', Ev' and Fermi levels (eV):
   electrons s->g, g->s, holes s->g, g->s */
static void tun_flux(const Mat *ms, const TPath *P, double Ecp, double Evp, double EFn, double EFp,
                     double EFm, double *F){
    double kT=ms->VT, h3=HPLANCK*HPLANCK*HPLANCK;
    double ce=4.*M_PI*ms->ms_e*M0*Q*Q*kT/h3, ch=4.*M_PI*ms->ms_h*M0*Q*Q*kT/h3;
    F[0]=ce*tun_int(P,Ecp,EFn,kT,0); F[1]=ce*tun_int(P,Ecp,EFm,kT,0);
    F[2]=ch*tun_int(P,Evp,EFp,kT,1); F[3]=ch*tun_int(P,Evp,EFm,kT,1);
}
/* one path: effective supply densities from the column, fluxes, and the
   linearised exchange: col_K (m^3/s per column node), tp_G (1/s) */
static void tun_path(Dev *d, int p, const double *Va){
    size_t k=d->tp_k[p]; int c=d->tp_c[p];
    const Mat *ms=MK(d,k); double VT=ms->VT, A=d->tp_A[p];
    TPath P; tun_stack(d,&P,d->tp_nl[p],d->tp_m+p*TUN_MAXL,d->tp_t+p*TUN_MAXL,d->phi[k],Va[c],d->cons[c].phi_m);
    int q0=d->tp_c0[p], nq=d->tp_cn[p];
    double Fs=P.D/(ms->eps_r*EPS0);                  /* surface field, + = pushes electrons to the wall */
    double neff[2], wsum[2];
    for(int hole=0;hole<2;hole++){
        double sg = hole ? -1. : 1., Fsc=sg*Fs;
        const double *car = hole ? d->p : d->n;
        double *K=d->col_K;
        for(int q=q0;q<q0+nq;q++) K[2*q+hole]=0.;
        if(Fsc>0.){                                  /* pressed against the wall: contact theorem */
            double P0=0., Fl=Fsc;
            for(int q=q0;q<q0+nq;q++){
                size_t kq=d->col_k[q];
                double Fr = (q<q0+nq-1 && d->col_h[q]>0.) ? sg*(d->phi[kq]-d->phi[d->col_k[q+1]])/d->col_h[q] : Fl;
                double Fm=0.5*(Fl+Fr); if(Fm<0.) Fm=0.;
                double w=Fm*d->col_l[q]/VT;          /* share of node q per unit density */
                K[2*q+hole]=w; P0+=w*car[kq];
                Fl=Fr;
            }
            size_t kl=d->col_k[q0+nq-1];
            K[2*(q0+nq-1)+hole]+=1.; P0+=car[kl];      /* kT n(z1) */
            neff[hole]=P0; wsum[hole]=1.;
        } else {                                     /* pushed away: local density at the interface */
            double qc = hole ? d->qcp[k] : d->qcn[k];
            neff[hole]=car[k]*exp(-qc+sg*(P.psis-d->phi[k])/VT);
            K[2*q0+hole]=neff[hole]/fmax(car[k],NFLOOR);
            wsum[hole]=1.;
        }
        if(!(neff[hole]>1e-100)) neff[hole]=1e-100;
        (void)wsum;
    }
    double EFn=-d->phi_n[k]-d->Cref, EFp=-d->phi_p[k]-d->Cref, EFm=-Va[c]-d->Cref;
    double Ecp=EFn-VT*log(neff[0]/ms->Nc), Evp=EFp+VT*log(neff[1]/ms->Nv);
    double F[4]; tun_flux(ms,&P,Ecp,Evp,EFn,EFp,EFm,F);
    /* K_q = A F_sg/n_eff * (d n_eff/d c_q) */
    double be=A*F[0]/neff[0], bh=A*F[2]/neff[1];
    for(int q=q0;q<q0+nq;q++){ d->col_K[2*q]*=be; d->col_K[2*q+1]*=bh; }
    d->tp_G[2*p]=A*F[1]; d->tp_G[2*p+1]=A*F[3];
}
static void tun_update(Dev *d, const double *Va){
    int np=d->n_tp; if(np<=0) return;
    #pragma omp parallel for schedule(dynamic,16) if(np>64)
    for(int p=0;p<np;p++) tun_path(d,p,Va);
}
static void tun_scatter(Dev *d, int hole){
    size_t N=NN(d);
    memset(d->tk,0,N*sizeof(double)); memset(d->tg,0,N*sizeof(double));
    for(int p=0;p<d->n_tp;p++){
        for(int q=d->tp_c0[p];q<d->tp_c0[p]+d->tp_cn[p];q++) d->tk[d->col_k[q]]+=d->col_K[2*q+hole];
        d->tg[d->tp_k[p]]+=d->tp_G[2*p+hole];
    }
}
/* escape rates (1/s) of path p at the current densities: electrons, holes */
static inline void tun_escape(const Dev *d, int p, double *se, double *sh){
    double a=0.,b=0.;
    for(int q=d->tp_c0[p];q<d->tp_c0[p]+d->tp_cn[p];q++){
        size_t kq=d->col_k[q]; a+=d->col_K[2*q]*d->n[kq]; b+=d->col_K[2*q+1]*d->p[kq];
    }
    *se=a; *sh=b;
}

/* Residual floor of the continuity solves, relative to sum|row terms|.
   Below ~1e-14 BiCGSTAB chases round-off, which floating regions (a MOS
   inversion layer fed only by generation) amplify into wrong solutions
   (verified: 2e-15 breaks the MOS-C test).  The price: smooth error modes of
   long quasi-neutral regions keep du ~ 1e-9..1e-8, i.e. a majority-carrier
   current noise ~sigma*VT*du/h -- the terminal-current floor reported below. */
static double CONT_ATOL = 2e-14;
static int solve_cont(Dev *d, int hole, double *dq, double lo, double hi, double *dmax){
    size_t N=NN(d); int its;
    tun_scatter(d,hole);
    cont_assemble(d,hole);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++) d->sol[k]=1.;
    double atol = CONT_ATOL;   /* residual floor relative to the row magnitudes */
    double rn=bicgstab(d,d->rhs,d->sol,d->max_iter_c>0?d->max_iter_c:1000,d->bicg_rtol,atol,&its);
    if(d->verbose>2) fprintf(stderr,"    bicgstab(%c) its=%d |r|=%.3e\n",hole?'p':'n',its,rn);
    int ok=isfinite(rn);
    if(ok) float_fix(d,hole);
    *dmax=cont_update(d,hole,dq,lo,hi);
    return ok;
}

/* Terminal currents: flux out of each contact's Dirichlet nodes into
   semiconductor nodes not belonging to the same contact (A, into device). */
static void terminal_currents(Dev *d, double *I, double *S){
    int nc=d->n_cons;
    for(int c=0;c<nc;c++){ I[c]=0.; S[c]=0.; d->cons[c].In=0.; }
    const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz;
    if(d->unstr) terminal_flux_g(d,I,S);
    else
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        int c=d->con_mask[k]; if(!c || d->kind[k]!=K_SEMI) continue;
        for(int dir=0;dir<6;dir++){
            size_t kb; double h,ar,Jn,Jp,Sn,Sp;
            if(!nbr(d,i,j,l,dir,&kb,&h,&ar) || d->kind[kb]!=K_SEMI || d->con_mask[kb]==c) continue;
            face_flux(d,k,kb,h,ar,&Jn,&Jp,&Sn,&Sp);
            I[c-1]+=Jn+Jp; S[c-1]+=Sn+Sp; d->cons[c-1].In+=Jn;
        }
    }
    /* gate tunnelling: electrons leaving to the gate = current into the
       device through the gate; holes leaving = current out of it */
    for(int p=0;p<d->n_tp;p++){
        int c=d->tp_c[p]; const double *G=d->tp_G+2*p; double se,sh;
        tun_escape(d,p,&se,&sh);
        double Ie=Q*(se-G[0]), Ih=-Q*(sh-G[1]);
        I[c]+=Ie+Ih; S[c]+=Q*(se+G[0]+sh+G[1]); d->cons[c].In+=Ie;
    }
}

/* Node current density (A/m^2) for display: mean of the adjacent face
   fluxes along each axis; recombination rate for display. */
static void upd_J(Dev *d){
    const int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz;
    if(d->unstr) upd_J_g(d);
    else
    #pragma omp parallel for collapse(2) schedule(static) if(NN(d)>PAR_MIN)
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        double *Jo[3]={d->Jx,d->Jy,d->Jz};
        if(d->kind[k]!=K_SEMI){ d->Jx[k]=d->Jy[k]=d->Jz[k]=0.; d->R_tot[k]=0.; continue; }
        for(int ax=0;ax<3;ax++){
            double s=0.; int cnt=0;
            size_t kb; double h,ar,Jn,Jp,Sn,Sp;
            if(nbr(d,i,j,l,2*ax,&kb,&h,&ar) && d->kind[kb]==K_SEMI){
                face_flux(d,kb,k,h,ar,&Jn,&Jp,&Sn,&Sp); s+=(Jn+Jp)/ar; cnt++; }
            if(nbr(d,i,j,l,2*ax+1,&kb,&h,&ar) && d->kind[kb]==K_SEMI){
                face_flux(d,k,kb,h,ar,&Jn,&Jp,&Sn,&Sp); s+=(Jn+Jp)/ar; cnt++; }
            Jo[ax][k]= cnt? s/cnt : 0.;
        }
        double R,a,b; recombK(d,k,&R,&a,&b);
        d->R_tot[k]=R;
    }
    /* gate tunnelling current density at the interface nodes (A/m^2, into
       the device through the gate) - kept in tg until the next solve */
    size_t N=NN(d);
    memset(d->tk,0,N*sizeof(double)); memset(d->tg,0,N*sizeof(double));
    for(int p=0;p<d->n_tp;p++){
        size_t k=d->tp_k[p]; const double *G=d->tp_G+2*p; double se,sh;
        tun_escape(d,p,&se,&sh);
        d->tg[k]+=Q*((se-G[0])-(sh-G[1])); d->tk[k]+=d->tp_A[p];
    }
    for(size_t k=0;k<N;k++) if(d->tk[k]>0.) d->tg[k]/=d->tk[k];
}

/* ===================================================== Gummel loop */
static void set_vrange(Dev *d, const double *Va){
    int any=0; double lo=0.,hi=0.;
    for(int c=0;c<d->n_cons;c++){
        if(d->cons[c].bc==2 && !d->cons[c].ntp) continue;   /* a tunnelling gate is a carrier terminal */
        if(!any){ lo=hi=Va[c]; any=1; }
        else { if(Va[c]<lo)lo=Va[c]; if(Va[c]>hi)hi=Va[c]; }
    }
    d->vlo=lo; d->vhi=hi;
}

static double newton_tol(const Dev *d, const double *Va){
    double s=1.;
    for(int c=0;c<d->n_cons;c++) if(fabs(Va[c])>s) s=fabs(Va[c]);
    return 2e-14+2e-15*s;
}

/* Anderson acceleration (Walker & Ni) of the full Gummel map
   x=(psi,phi_n,phi_p) -> g(x), depth d->aa_m (runtime cap AA_DEPTH).
   Mixing phi_n,p as well as psi matters: the carrier densities and the
   field-dependent mobility of iteration k+1 are built from phi_n,p, so a
   psi-only mixing sees a map with hidden state (erratic in saturated
   channels).  Least squares by regularised normal equations on the cached
   Gram matrix of the residual differences.
   phi_n (phi_p) enters the least squares and the mixing only where that
   carrier matters (> 1e-6 of n+p+|Nd-Na|, the same test as the convergence
   check): elsewhere the quasi-Fermi potential of a vanishing carrier is
   numerically meaningless, can move by volts per iteration, and would
   swamp the residual norm; those nodes take the plain Gummel update. */
static inline int aa_sig(const Dev *d, size_t k, int hole){
    if(d->kind[k]!=K_SEMI || d->con_mask[k]) return 0;
    double c = hole ? d->p[k] : d->n[k];
    return c > 1e-6*(d->n[k]+d->p[k]+fabs(d->Nd[k]-d->Na[k]));
}
static int AA_DEPTH = AA_MAX;    /* env SEMISIM_AAM (0 = plain Gummel) */
static void aa_reset(Dev *d){ d->aa_n=0; d->aa_pos=0; d->aa_have=0; }
static double dot_fd(const float *a, const double *b, size_t L){
    double s=0.;
    #pragma omp parallel for reduction(+:s) schedule(static) if(L>PAR_MIN)
    for(size_t k=0;k<L;k++) s+=(double)a[k]*b[k];
    return s;
}
static double dot_ff(const float *a, const float *b, size_t L){
    double s=0.;
    #pragma omp parallel for reduction(+:s) schedule(static) if(L>PAR_MIN)
    for(size_t k=0;k<L;k++) s+=(double)a[k]*(double)b[k];
    return s;
}
/* x: the 3 input blocks; g: the 3 image blocks, overwritten by the mixed iterate */
static void aa_step(Dev *d, double *const x[3], double *const g[3]){
    size_t N=NN(d), L=3*N;
    int M = d->aa_m < AA_DEPTH ? d->aa_m : AA_DEPTH;
    double *f=d->aa_f, *gp=d->aa_g;
    if(M>0 && d->aa_have){
        int s=d->aa_pos; float *dF=d->aa_dF[s], *dG=d->aa_dG[s];
        for(int b=0;b<3;b++){
            const double *xb=x[b], *gb=g[b];
            #pragma omp parallel for schedule(static) if(N>PAR_MIN)
            for(size_t k=0;k<N;k++){ size_t q=b*N+k;
                double fk = (b==0 || aa_sig(d,k,b-1)) ? gb[k]-xb[k] : 0.;
                dF[q]=(float)(fk-f[q]); dG[q]=(float)(gb[k]-gp[q]); }
        }
        if(d->aa_n<M) d->aa_n++;
        for(int i=0;i<d->aa_n;i++){                 /* refresh Gram row/col s */
            double v=dot_ff(d->aa_dF[i],dF,L); d->aa_gram[i][s]=d->aa_gram[s][i]=v; }
        d->aa_pos=(s+1)%M;
    }
    for(int b=0;b<3;b++){
        const double *xb=x[b], *gb=g[b];
        #pragma omp parallel for schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){ size_t q=b*N+k;
            f[q] = (b==0 || aa_sig(d,k,b-1)) ? gb[k]-xb[k] : 0.; gp[q]=gb[k]; }
    }
    d->aa_have=1;
    int m=d->aa_n; if(m==0) return;                 /* plain Picard: x = g */
    double A[AA_MAX][AA_MAX], bb[AA_MAX], gam[AA_MAX];
    for(int i=0;i<m;i++){ bb[i]=dot_fd(d->aa_dF[i],f,L); for(int j=0;j<m;j++) A[i][j]=d->aa_gram[i][j]; }
    double tr=0.; for(int i=0;i<m;i++) tr+=A[i][i];
    for(int i=0;i<m;i++) A[i][i]+=1e-9*tr+1e-300;   /* float differences: ~1e-7 relative */
    for(int c=0;c<m;c++){                            /* Gaussian elimination */
        int pv=c; for(int r=c+1;r<m;r++) if(fabs(A[r][c])>fabs(A[pv][c])) pv=r;
        if(pv!=c){ for(int q=0;q<m;q++){ double t=A[c][q];A[c][q]=A[pv][q];A[pv][q]=t; } double t=bb[c];bb[c]=bb[pv];bb[pv]=t; }
        for(int r=c+1;r<m;r++){ double fct=A[r][c]/A[c][c]; for(int q=c;q<m;q++) A[r][q]-=fct*A[c][q]; bb[r]-=fct*bb[c]; }
    }
    for(int r=m-1;r>=0;r--){ double t=bb[r]; for(int q=r+1;q<m;q++) t-=A[r][q]*gam[q]; gam[r]=t/A[r][r]; }
    for(int i=0;i<m;i++) if(!(gam[i]==gam[i])) return;   /* NaN: keep Picard */
    for(int b=0;b<3;b++){
        double *gb=g[b];
        #pragma omp parallel for schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){
            if(d->con_mask[k]) continue;
            if(b>0 && !aa_sig(d,k,b-1)) continue;  /* plain Gummel update */
            double s=gb[k]; size_t q=b*N+k;
            for(int i=0;i<m;i++) s-=gam[i]*(double)d->aa_dG[i][q];
            gb[k]=s;
        }
    }
}

/* Quasi-Fermi potentials are clamped to the span of the carrier-terminal
   voltages: in steady state without external generation, drift-diffusion
   with SRH/Auger/radiative recombination obeys a maximum principle
   (min V_c <= phi_n,phi_p <= max V_c).  The exact bound also pins floating
   regions fed only by thermal generation (a MOS inversion layer) whose
   continuity block is numerically near-singular. */
static double QF_MARGIN = 0.;
static void clamp_qf(Dev *d);
/* GUI hooks, touched from another thread while a solve runs: a cooperative
   stop request (checked once per Gummel iteration; the solver then falls
   back to the last converged bias point) and a progress record:
   0 continuation fraction t of the running ramp, 1 Gummel iteration,
   2 residual (V), 3 accepted continuation steps, 4 phase (0 idle,
   1 equilibrium, 2 gate ramp, 3 bias ramp), 5 total Gummel iterations. */
static volatile int g_stop = 0;
static volatile double g_prog[8];
/* One Gummel solve at fixed terminal voltages.
   Iteration k, state x_k=(psi_k,phi_n,k,phi_p,k):
     carriers + field-dependent mobility from x_k -> continuity at psi_k
     (currents evaluated here, exactly conservative) -> phi_n*,phi_p* ->
     Poisson -> psi* ;  g(x_k)=(psi*,phi_n*,phi_p*) -> Anderson mix -> x_{k+1}.
   The accepted state is psi_k with its own continuity solution. */
static int gummel(Dev *d, const double *Va, int maxit, double gtol, int *iters, double *res){
    size_t N=NN(d); int nc=d->n_cons;
    double lo=d->vlo-QF_MARGIN, hi=d->vhi+QF_MARGIN, nt0=newton_tol(d,Va);
    double I[NC_MAX],S[NC_MAX],Iold[NC_MAX]; int haveI=0;
    double *xk=d->wk[6], *dqn=d->Jx, *dqp=d->Jy;
    double *qn0=d->Jz, *qp0=d->R_tot;       /* phi_n,p entering the iteration (display arrays are free here) */
    double itol = d->tol_c>0. ? d->tol_c : 1e-6;
    apply_bc(d,Va);
    aa_reset(d);
    /* start from a Poisson-consistent potential for the predicted phi_n,p */
    poisson_newton(d,Va,fmax(nt0,1e-6),100);
    int it, conv=0, it_mark=0; double dpsi=1.,dq=1., best=1e300, prev=1., best_mark=1e300;
    for(it=1;it<=maxit;it++){
        if(g_stop) break;
        double t0=wtime();
        clamp_qf(d);
        memcpy(qn0,d->phi_n,N*sizeof(double)); memcpy(qp0,d->phi_p,N*sizeof(double));
        rebuild_np(d);
        upd_mob(d);
        tun_update(d,Va);
        double t1=wtime(); d->t_mob+=t1-t0;
        double dn=0.,dp=0.;
        if(!solve_cont(d,0,dqn,lo,hi,&dn)) break;
        if(!solve_cont(d,1,dqp,lo,hi,&dp)) break;
        double t2=wtime(); d->t_cont+=t2-t1;
        dq=fmax(dn,dp);
        terminal_currents(d,I,S);
        d->t_cur+=wtime()-t2;
        int icv=haveI;
        if(haveI) for(int c=0;c<nc;c++)
            if(fabs(I[c]-Iold[c])>fmax(itol*fabs(I[c]),3e-14*S[c])) icv=0;
        memcpy(Iold,I,sizeof(double)*nc); haveI=1;
        /* Poisson from the new quasi-Fermi potentials */
        memcpy(xk,d->phi,N*sizeof(double));
        #pragma omp parallel for schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){                 /* neutral-region predictor */
            if(d->kind[k]!=K_SEMI || d->con_mask[k]) continue;
            double n=d->n[k],p=d->p[k],D=fabs(d->Nd[k]-d->Na[k]);
            double den=fmax(n+p,D); if(!(den>0.)) continue;
            d->phi[k]+=(n*dqn[k]+p*dqp[k])/den;
        }
        double ntol=fmax(nt0,d->ntol_fac*fmin(prev,1.));
        double t3=wtime();
        poisson_newton(d,Va,ntol,100);
        d->t_pois+=wtime()-t3;
        dpsi=0.;
        #pragma omp parallel for reduction(max:dpsi) schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){ double a=fabs(d->phi[k]-xk[k]); if(a>dpsi) dpsi=a; }
        if(d->verbose>1)
            fprintf(stderr,"  gummel %3d dpsi=%.3e dq=%.3e I0=%.6e aa=%d\n",it,dpsi,dq,nc?I[0]:0.,d->aa_n);
        if(dpsi!=dpsi || dq!=dq) break;
        double r=fmax(dpsi,dq);
        g_prog[1]=it; g_prog[2]=r; g_prog[5]=(double)(d->gummel_its+it);
        if(dq<gtol && dpsi<gtol && icv){           /* accept psi_k + its carriers */
            memcpy(d->phi,xk,N*sizeof(double)); conv=1; break;
        }
        prev=r;
        /* slow but steady convergence is fine (high-current Gummel); stop
           only on stagnation: best residual not halved within 60 iterations.
           A stagnating iterate is still accepted if it is within 1e3*gtol and
           every terminal current is stable to tol_c (or below its round-off
           floor): that is the noise floor of a floating region (the n- base
           of a blocking thyristor, fed only by generation), whose majority
           quasi-Fermi level double precision cannot pin any tighter. */
        if(r<0.5*best_mark){ best_mark=r; it_mark=it; }
        if(it-it_mark>20 && r<=1e3*gtol && icv){
            memcpy(d->phi,xk,N*sizeof(double)); conv=1;
            if(d->verbose) fprintf(stderr,"semisim: accepted at noise floor (res %.2e)\n",r);
            break;
        }
        if(it-it_mark>60) break;
        if(r<best) best=r;
        else if(r>30.*best) aa_reset(d);          /* AA went astray: restart */
        double t4=wtime();
        double *xx[3]={xk,qn0,qp0}, *gg[3]={d->phi,d->phi_n,d->phi_p};
        aa_step(d,xx,gg);
        d->t_aa+=wtime()-t4;
    }
    if(it>maxit) it=maxit;
    d->gummel_its+=it; *iters=it; *res=fmax(dpsi,dq);
    return conv;
}

/* ================================================= continuation */
static void save_state(Dev *d, double **h){
    size_t N=NN(d);
    memcpy(h[0],d->phi,N*sizeof(double));
    memcpy(h[1],d->phi_n,N*sizeof(double));
    memcpy(h[2],d->phi_p,N*sizeof(double));
}
static void load_state(Dev *d, double **h){
    size_t N=NN(d);
    memcpy(d->phi,h[0],N*sizeof(double));
    memcpy(d->phi_n,h[1],N*sizeof(double));
    memcpy(d->phi_p,h[2],N*sizeof(double));
}
static void clamp_qf(Dev *d){
    size_t N=NN(d); double lo=d->vlo-QF_MARGIN, hi=d->vhi+QF_MARGIN;
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        if(d->phi_n[k]<lo) d->phi_n[k]=lo; else if(d->phi_n[k]>hi) d->phi_n[k]=hi;
        if(d->phi_p[k]<lo) d->phi_p[k]=lo; else if(d->phi_p[k]>hi) d->phi_p[k]=hi;
    }
}

/* First-step predictor: the quasi-neutral region attached to each ohmic
   contact moves rigidly with that contact's voltage change. */
static void region_shift(Dev *d, const double *dV){
    if(d->unstr){ region_shift_g(d,dV); return; }
    size_t N=NN(d);
    unsigned char *mark=calloc(N,1); if(!mark) return;
    size_t *q=(size_t*)d->wk[0];
    for(int c=0;c<d->n_cons;c++){
        if(d->cons[c].bc!=0 || fabs(dV[c])<1e-15) continue;
        size_t head=0,tail=0; int type=0;
        for(size_t k=0;k<N;k++) if(d->con_mask[k]==c+1 && d->kind[k]==K_SEMI && !mark[k]){
            mark[k]=(unsigned char)(c+1); q[tail++]=k;
            if(!type){ double D=d->Nd[k]-d->Na[k]; type = D>0.?1:(D<0.?-1:0); }
        }
        if(!type) continue;
        while(head<tail){
            size_t k=q[head++]; int i,j,l; unidx(d,k,&i,&j,&l);
            for(int dir=0;dir<6;dir++){
                size_t kb; double h,ar;
                if(!nbr(d,i,j,l,dir,&kb,&h,&ar) || mark[kb] || d->kind[kb]!=K_SEMI) continue;
                if(d->con_mask[kb] && d->con_mask[kb]!=c+1) continue;
                double D=d->Nd[kb]-d->Na[kb]; if(D*type<=0.) continue;
                double maj = type>0 ? d->n[kb] : d->p[kb];
                if(maj<0.5*fabs(D)) continue;
                mark[kb]=(unsigned char)(c+1); q[tail++]=kb;
            }
        }
        for(size_t t=0;t<tail;t++){
            size_t k=q[t]; if(d->con_mask[k]) continue;
            d->phi[k]+=dV[c]; d->phi_n[k]+=dV[c]; d->phi_p[k]+=dV[c];
        }
    }
    free(mark);
}

static int equilibrium(Dev *d){
    double Va[NC_MAX];
    for(int c=0;c<d->n_cons;c++) Va[c] = d->cons[c].bc==2 ? d->cons[c].V : 0.;
    init_eq(d);
    mob0_init(d);
    set_vrange(d,Va);
    apply_bc(d,Va);
    g_prog[4]=1.;
    int ok=poisson_newton(d,Va,newton_tol(d,Va),300);
    rebuild_np(d); upd_mob(d);
    for(int c=0;c<d->n_cons;c++) d->cons[c].Vprev=Va[c];
    d->have_solution=1;
    if(d->verbose) fprintf(stderr,"semisim: equilibrium %s (newton its %ld)\n",ok?"ok":"NOT converged",d->newton_its);
    return ok;
}

static double gummel_tol(const Dev *d){
    double g=d->tol_p; if(!(g>0.)) g=1e-9;
    if(g>1e-4) g=1e-4; if(g<1e-13) g=1e-13;
    return g;
}

/* Adaptive continuation from the voltages of the last solution (Vprev) to
   the requested ones (V), all contacts ramped together. */
static int cont_path(Dev *d){
    int nc=d->n_cons; double Vf[NC_MAX],Vt[NC_MAX],Va[NC_MAX],dV[NC_MAX];
    double span=0.;
    for(int c=0;c<nc;c++){ Vf[c]=d->cons[c].Vprev; Vt[c]=d->cons[c].V;
                           if(fabs(Vt[c]-Vf[c])>span) span=fabs(Vt[c]-Vf[c]); }
    int maxg = d->max_gummel>0 ? d->max_gummel : 100;
    double gtol=gummel_tol(d);
    size_t N=NN(d);
    mob0_init(d); upd_mob(d);
    save_state(d,d->h1);
    double t=0., tp=-1., dt = span>0. ? fmin(1.,0.1/span) : 1.;
    int ok_all=0, its=0, steps=0, fails=0; double res=0.;
    g_prog[0]=0.; g_prog[3]=0.;
    for(int guard=0;guard<5000;guard++){
        if(g_stop) break;
        double tn=fmin(1.,t+dt);
        for(int c=0;c<nc;c++) Va[c]=Vf[c]+tn*(Vt[c]-Vf[c]);
        set_vrange(d,Va);
        if(tp>=0.){
            double r=(tn-t)/(t-tp);
            #pragma omp parallel for schedule(static) if(N>PAR_MIN)
            for(size_t k=0;k<N;k++){
                for(int a=0;a<3;a++){
                    double x1=d->h1[a][k], x2=d->h2[a][k];
                    double *dst = a==0?d->phi:(a==1?d->phi_n:d->phi_p);
                    dst[k]=x1+r*(x1-x2);
                }
            }
        } else {
            load_state(d,d->h1);
            if(span>0.){
                for(int c=0;c<nc;c++) dV[c]=(tn-t)*(Vt[c]-Vf[c]);
                region_shift(d,dV);
            }
        }
        clamp_qf(d); apply_bc(d,Va); rebuild_np(d);
        int ok=gummel(d,Va,maxg,gtol,&its,&res);
        d->conv_iter=its; d->conv_res=res;
        if(d->verbose)
            fprintf(stderr,"semisim: step t=%.4f (V0=%.4f) %s its=%d res=%.2e\n",
                    tn,nc?Va[0]:0.,ok?"ok":"FAIL",its,res);
        if(ok){
            double *tmp[3]={d->h2[0],d->h2[1],d->h2[2]};
            for(int a=0;a<3;a++){ d->h2[a]=d->h1[a]; d->h1[a]=tmp[a]; }
            save_state(d,d->h1);
            tp=t; t=tn; steps++;
            g_prog[0]=t; g_prog[3]=steps;
            for(int c=0;c<nc;c++) d->cons[c].Vprev=Va[c];
            if(t>=1.){ ok_all=1; break; }
            /* grow on easy steps; never shrink after a success (with a
               linearly converging Gummel map the iteration count is set by
               the contraction rate, not by the step length) */
            dt*= its<=8 ? 2.0 : (its<=25 ? 1.5 : 1.2);
        } else {
            fails++;
            if(span==0. || dt*span<2e-5){ break; }
            dt*=0.3;
        }
    }
    if(!ok_all && (span>0. || g_stop)){
        /* keep the last converged state (partial bias) */
        load_state(d,d->h1);
        for(int c=0;c<nc;c++) Va[c]=d->cons[c].Vprev;
        set_vrange(d,Va); apply_bc(d,Va); rebuild_np(d); upd_mob(d);
        if(d->verbose) fprintf(stderr,"semisim: continuation stopped at t=%.4f\n",t);
    }
    return ok_all;
}

/* Bias path: when MOS gates and carrier terminals both change, ramp the
   gates first (carrier terminals held), then the terminals - the order of a
   real measurement.  Ramping both together drags e.g. an IGBT through
   channel formation and the onset of high injection at the same time. */
static int run_continuation(Dev *d){
    int nc=d->n_cons, gch=0, cch=0, ok_all;
    for(int c=0;c<nc;c++){
        if(fabs(d->cons[c].V-d->cons[c].Vprev)<=1e-12) continue;
        if(d->cons[c].bc==2) gch=1; else cch=1;
    }
    if(gch && cch){
        double Vt[NC_MAX];
        for(int c=0;c<nc;c++){ Vt[c]=d->cons[c].V; if(d->cons[c].bc!=2) d->cons[c].V=d->cons[c].Vprev; }
        if(d->verbose) fprintf(stderr,"semisim: gate ramp first\n");
        g_prog[4]=2.;
        ok_all=cont_path(d);
        for(int c=0;c<nc;c++) d->cons[c].V=Vt[c];
        g_prog[4]=3.;
        if(ok_all) ok_all=cont_path(d);
    } else { g_prog[4]= gch ? 2. : 3.; ok_all=cont_path(d); }
    double I[NC_MAX],S[NC_MAX];
    terminal_currents(d,I,S);
    /* Terminal currents and their resolution floors.  The flux out of a
       contact is a small difference of large drift and diffusion terms when
       the contact region is heavily doped: its round-off floor is ~3e-14 of
       the flux magnitudes S (the continuity solves resolve relative densities
       to ~1e-14).  By Kirchhoff the same current is also minus the sum of the
       other terminals, whose floors may be far lower (a p+ anode against a
       lightly doped cathode): report whichever estimate is better resolved.
       If the observed mismatch exceeds the floors, scale them up (honest
       error bars for a partially converged state). */
    /* gates (tunnelling current, 0 without it) take part in the sum but
       always report their own, directly evaluated current */
    double ks=0., F=0., fl[NC_MAX]; int ncar=0;
    for(int c=0;c<nc;c++){ ks+=I[c]; fl[c]=3e-14*S[c]; F+=fl[c]; if(d->cons[c].bc!=2) ncar++; }
    double scl = F>0. ? fmax(1.,2.*fabs(ks)/F) : 1.;
    for(int c=0;c<nc;c++){
        double Iv=I[c], fv=fl[c], fo=F-fl[c];
        if(d->cons[c].bc!=2 && ncar>1 && fo<fv){ Iv=-(ks-I[c]); fv=fo; }
        d->cons[c].I=Iv; d->cons[c].Iflr=fv*scl; d->cons[c].Iraw=I[c];
    }
    return ok_all;
}

/* ============================================ unstructured (graph) transport
   The same finite-volume equations as the tensor code, with the neighbours
   and the box-method geometry taken from the CSR graph (d->unstr = 1):
   edge length gh, dual face area ga, control volume gvol.  Matrix
   off-diagonals are stored per directed edge (gval), the diagonal in A.c. */

static void poisson_core_g(Dev *d, const double *phi, int build, double *F){
    size_t N=NN(d); double *dgv=d->A.c;
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        int e0=d->grp[k], e1=d->grp[k+1];
        if(is_dir(d,k)){
            F[k]=0.;
            if(build){ dgv[k]=1.; for(int e=e0;e<e1;e++) d->gval[e]=0.; d->sc[k]=1.; }
            continue;
        }
        double ek=d->eps[k], pk=phi[k], f=0., dg=0.;
        for(int e=e0;e<e1;e++){
            size_t kb=d->gcol[e]; double eb=d->eps[kb], cf=(2.*ek*eb/(ek+eb))*d->ga[e]/d->gh[e];
            f+=cf*(phi[kb]-pk); dg+=cf;
            if(build) d->gval[e]= is_dir(d,kb) ? 0. : -cf;
        }
        double vol=d->gvol[k];
        if(d->kind[k]==K_SEMI){
            const Mat *m=&d->mats[d->mat_idx[k]];
            double n=nK(d,k,pk,d->phi_n[k]), p=pK(d,k,pk,d->phi_p[k]);
            f +=vol*Q*(p-n+d->Nd[k]-d->Na[k]+d->Nf[k]);
            dg+=vol*Q*(n+p)/m->VT;
        } else f+=vol*Q*d->Nf[k];
        F[k]=f;
        if(build){ dgv[k]=dg; d->sc[k]=dg; }
    }
}

static inline void face_flux_g(const Dev *d, size_t ka, size_t kb, int ue, double h, double ar,
                               double *Jn, double *Jp, double *Sn, double *Sp){
    const Mat *ma=MK(d,ka);
    double VT=ma->VT, dps=(d->phi[kb]-d->phi[ka])/VT;
    double Dn=dps+dlnc(d,ka,kb), Dp=dps-dlnv(d,ka,kb);
    double gn=Q*ar*(double)d->gmu[0][ue]*VT/h, gp=Q*ar*(double)d->gmu[1][ue]*VT/h;
    double a1=gn*d->n[kb]*B(Dn), a2=gn*d->n[ka]*B(-Dn);
    double b1=gp*d->p[ka]*B(Dp), b2=gp*d->p[kb]*B(-Dp);
    *Jn=a1-a2; *Jp=b1-b2; *Sn=a1+a2; *Sp=b1+b2;
}

/* edge mobilities (Lombardi: the field normal to the edge from the node
   gradients, |grad psi|^2 - (grad psi . t)^2) */
static void upd_mob_ref_g(Dev *d, const double *qn_, const double *qp_){
    size_t N=NN(d);
    const int sat=(d->models&4)!=0, lb=(d->models&16)!=0;
    if(lb){
        #pragma omp parallel for schedule(static) if(N>PAR_MIN)
        for(size_t k=0;k<N;k++){
            if(d->kind[k]==K_SEMI && d->qd[k]<0.5) g_grad(d,d->phi,k,d->ggrad+3*k);
            else { d->ggrad[3*k]=d->ggrad[3*k+1]=d->ggrad[3*k+2]=0.; }
        }
    }
    #pragma omp parallel for schedule(dynamic,256) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        for(int e=d->grp[k];e<d->grp[k+1];e++){
            size_t kp=d->gcol[e]; if(kp<k) continue;          /* each undirected edge once */
            int ue=d->gue[e];
            if(d->kind[k]!=K_SEMI || d->kind[kp]!=K_SEMI){ d->gmu[0][ue]=d->gmu[1][ue]=0.; continue; }
            const Mat *m=MK(d,k); double h=d->gh[e];
            double mn=0.5*(d->mu_n[k]+d->mu_n[kp]), mp=0.5*(d->mu_p[k]+d->mu_p[kp]);
            if(lb && m->lb_on){
                double Dw=exp(-0.5*(d->qd[k]+d->qd[kp])/LB_LCRIT);
                if(Dw>1e-8){
                    double t[3]; for(int q=0;q<3;q++) t[q]=(d->gpos[3*kp+q]-d->gpos[3*k+q])/h;
                    const double *ga_=d->ggrad+3*k, *gb=d->ggrad+3*kp;
                    double pa=ga_[0]*t[0]+ga_[1]*t[1]+ga_[2]*t[2], pb=gb[0]*t[0]+gb[1]*t[1]+gb[2]*t[2];
                    double Fa=ga_[0]*ga_[0]+ga_[1]*ga_[1]+ga_[2]*ga_[2]-pa*pa;
                    double Fb=gb[0]*gb[0]+gb[1]*gb[1]+gb[2]*gb[2]-pb*pb;
                    double F=sqrt(fmax(0.5*(Fa+Fb),0.));
                    double Nt=0.5*(d->Nd[k]+d->Na[k]+d->Nd[kp]+d->Na[kp]);
                    mn=lombardi(m,mn,F,Nt,Dw,0); mp=lombardi(m,mp,F,Nt,Dw,1);
                }
            }
            if(sat){
                double dps=d->phi[kp]-d->phi[k];
                double Fn=heat_field(qn_[kp]-qn_[k],dps,h), Fp=heat_field(qp_[kp]-qp_[k],dps,h);
                mn=mob_hf(mn,Fn,m->vsat_n,m->beta_n); mp=mob_hf(mp,Fp,m->vsat_p,m->beta_p);
            }
            d->gmu[0][ue]=(mu_t)mn; d->gmu[1][ue]=(mu_t)mp;
        }
    }
}

static void mob_display_g(Dev *d){
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        if(d->kind[k]!=K_SEMI) continue;
        double sn=0.,sp=0.; int c=0;
        for(int e=d->grp[k];e<d->grp[k+1];e++){
            if(d->kind[d->gcol[e]]!=K_SEMI) continue;
            sn+=d->gmu[0][d->gue[e]]; sp+=d->gmu[1][d->gue[e]]; c++;
        }
        if(c){ d->wk[0][k]=sn/c; d->wk[1][k]=sp/c; }
        else { d->wk[0][k]=d->mu_n[k]; d->wk[1][k]=d->mu_p[k]; }
    }
    for(size_t k=0;k<N;k++) if(d->kind[k]==K_SEMI){ d->mu_n[k]=d->wk[0][k]; d->mu_p[k]=d->wk[1][k]; }
}

static void cont_assemble_g(Dev *d, int hole){
    size_t N=NN(d); double *dgv=d->A.c;
    const double *car = hole? d->p : d->n;
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        int e0=d->grp[k], e1=d->grp[k+1];
        if(d->kind[k]!=K_SEMI || d->con_mask[k]){
            dgv[k]=1.; for(int e=e0;e<e1;e++) d->gval[e]=0.;
            d->rhs[k]=(d->kind[k]==K_SEMI)?1.:0.; d->sc[k]=(d->kind[k]==K_SEMI)?car[k]:0.;
            d->dsc[k]=0.;
            continue;
        }
        const Mat *ma=MK(d,k);
        double VT=ma->VT, sk=fmax(car[k],NFLOOR), dg=0., b=0.;
        for(int e=e0;e<e1;e++){
            size_t kb=d->gcol[e]; d->gval[e]=0.;
            if(d->kind[kb]!=K_SEMI) continue;
            double D=(d->phi[kb]-d->phi[k])/VT + (hole ? -dlnv(d,k,kb) : dlnc(d,k,kb));
            double g=d->ga[e]*(double)d->gmu[hole][d->gue[e]]*VT/d->gh[e];
            double self = hole ? g*B(D)  : g*B(-D);
            double nb   = hole ? g*B(-D) : g*B(D);
            dg+=self;
            if(d->con_mask[kb]) b+=nb*car[kb];
            else d->gval[e]=-nb*fmax(car[kb],NFLOOR);
        }
        double R,dRn,dRp; recombK(d,k,&R,&dRn,&dRp);
        double dR = hole?dRp:dRn; if(dR<0.) dR=0.;
        double vol=d->gvol[k];
        dg+=vol*dR; b+=vol*(dR*car[k]-R);
        dg+=d->tk[k]; b+=d->tg[k];
        double w=1./(dg*sk);
        dgv[k]=1.;
        for(int e=e0;e<e1;e++) d->gval[e]*=w;
        d->rhs[k]=b*w; d->sc[k]=sk; d->dsc[k]=dg*sk;
    }
    d->mg_w=d->dsc; d->mg_qsign= hole ? -1 : 1;
}

/* floating majority regions (see float_fix) on the graph */
static int float_fix_g(Dev *d, int hole){
    size_t N=NN(d);
    int *lab=(int*)d->wk[0]; size_t *q=(size_t*)d->wk[1];
    const double *oth = hole ? d->n : d->p;
    double *u=d->sol, *s=d->sc;
    #define MAJ(k) (d->kind[k]==K_SEMI && !d->con_mask[k] && s[k]*u[k]>=oth[k] \
                    && s[k]*u[k]>=0.5*fabs(d->Nd[k]-d->Na[k]) && s[k]*u[k]>10.*niK(d,k))
    for(size_t k=0;k<N;k++) lab[k]=0;
    int nfix=0, blk=0;
    for(size_t k0=0;k0<N;k0++){
        if(lab[k0] || !MAJ(k0)) continue;
        blk++; size_t head=0,tail=0; q[tail++]=k0; lab[k0]=blk; int anch=0;
        while(head<tail){
            size_t k=q[head++];
            for(int e=d->grp[k];e<d->grp[k+1];e++){
                size_t kb=d->gcol[e];
                if(d->kind[kb]!=K_SEMI) continue;
                if(d->con_mask[kb]){ anch=1; continue; }
                if(!lab[kb] && MAJ(kb)){ lab[kb]=blk; q[tail++]=kb; }
            }
        }
        if(anch) continue;
        double Aout=0., Bin=0.;
        for(size_t t=0;t<tail;t++){
            size_t k=q[t]; const Mat *ma=MK(d,k); double VT=ma->VT;
            for(int e=d->grp[k];e<d->grp[k+1];e++){
                size_t kb=d->gcol[e];
                if(d->kind[kb]!=K_SEMI || lab[kb]==blk) continue;
                double D=(d->phi[kb]-d->phi[k])/VT + (hole ? -dlnv(d,k,kb) : dlnc(d,k,kb));
                double g=d->ga[e]*(double)d->gmu[hole][d->gue[e]]*VT/d->gh[e];
                double self = hole ? g*B(D) : g*B(-D), nb = hole ? g*B(-D) : g*B(D);
                double cb = d->con_mask[kb] ? (hole?d->p[kb]:d->n[kb]) : s[kb]*u[kb];
                Aout+=self*s[k]*u[k]; Bin+=nb*cb;
            }
        }
        double sl=0.;
        for(int it=0;it<60;it++){
            double lam=exp(sl), G=lam*Aout-Bin, dG=lam*Aout;
            for(size_t t=0;t<tail;t++){
                size_t k=q[t];
                double c=lam*s[k]*u[k], R,dRn,dRp;
                recomb(MK(d,k),niK(d,k),tauf(d,k), hole?oth[k]:c, hole?c:oth[k], &R,&dRn,&dRp);
                double vol=d->gvol[k];
                G+=vol*R; dG+=vol*(hole?dRp:dRn)*c;
                G+=d->tk[k]*c-d->tg[k]; dG+=d->tk[k]*c;
            }
            if(!(dG>0.)) break;
            double ds=-G/dG; if(ds>2.) ds=2.; else if(ds<-2.) ds=-2.;
            sl+=ds; if(fabs(ds)<1e-13) break;
        }
        if(sl==sl && fabs(sl)>1e-15 && fabs(sl)<200.){
            double lam=exp(sl);
            for(size_t t=0;t<tail;t++) u[q[t]]*=lam;
            nfix++;
        }
    }
    #undef MAJ
    return nfix;
}

static void terminal_flux_g(Dev *d, double *I, double *S){
    size_t N=NN(d);
    for(size_t k=0;k<N;k++){
        int c=d->con_mask[k]; if(!c || d->kind[k]!=K_SEMI) continue;
        for(int e=d->grp[k];e<d->grp[k+1];e++){
            size_t kb=d->gcol[e]; double Jn,Jp,Sn,Sp;
            if(d->kind[kb]!=K_SEMI || d->con_mask[kb]==c) continue;
            face_flux_g(d,k,kb,d->gue[e],d->gh[e],d->ga[e],&Jn,&Jp,&Sn,&Sp);
            I[c-1]+=Jn+Jp; S[c-1]+=Sn+Sp; d->cons[c-1].In+=Jn;
        }
    }
}

/* node current-density vector from the edge fluxes by least squares,
   J = (sum_e A_e t t^T)^-1 sum_e A_e t j_e; which: 0 total, 1 n, 2 p */
static void g_node_J(const Dev *d, size_t k, int which, double *J){
    double M[6]={0,0,0,0,0,0}, r[3]={0,0,0};
    J[0]=J[1]=J[2]=0.;
    if(d->kind[k]!=K_SEMI) return;
    for(int e=d->grp[k];e<d->grp[k+1];e++){
        size_t kb=d->gcol[e]; if(d->kind[kb]!=K_SEMI) continue;
        double h=d->gh[e], ar=d->ga[e], t[3], Jn,Jp,Sn,Sp;
        for(int q=0;q<3;q++) t[q]=(d->gpos[3*kb+q]-d->gpos[3*k+q])/h;
        face_flux_g(d,k,kb,d->gue[e],h,ar,&Jn,&Jp,&Sn,&Sp);
        double j = (which==1 ? Jn : (which==2 ? Jp : Jn+Jp))/(ar>0.?ar:1e-300);
        double w=ar;
        M[0]+=w*t[0]*t[0]; M[1]+=w*t[0]*t[1]; M[2]+=w*t[0]*t[2];
        M[3]+=w*t[1]*t[1]; M[4]+=w*t[1]*t[2]; M[5]+=w*t[2]*t[2];
        r[0]+=w*t[0]*j; r[1]+=w*t[1]*j; r[2]+=w*t[2]*j;
    }
    double tr=M[0]+M[3]+M[5], eps=1e-6*tr+1e-300;
    M[0]+=eps; M[3]+=eps; M[5]+=eps;
    double a=M[0],b=M[1],c=M[2],dd=M[3],e=M[4],f=M[5];
    double A=dd*f-e*e, Bq=-(b*f-c*e), Cc=b*e-c*dd, D=a*f-c*c, E=-(a*e-b*c), F=a*dd-b*b;
    double det=a*A+b*Bq+c*Cc; if(!(fabs(det)>0.)) return;
    J[0]=(A*r[0]+Bq*r[1]+Cc*r[2])/det;
    J[1]=(Bq*r[0]+D*r[1]+E*r[2])/det;
    J[2]=(Cc*r[0]+E*r[1]+F*r[2])/det;
}
static void upd_J_g(Dev *d){
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++){
        if(d->kind[k]!=K_SEMI){ d->Jx[k]=d->Jy[k]=d->Jz[k]=0.; d->R_tot[k]=0.; continue; }
        double J[3]; g_node_J(d,k,0,J);
        d->Jx[k]=J[0]; d->Jy[k]=J[1]; d->Jz[k]=J[2];
        double R,a,b; recombK(d,k,&R,&a,&b);
        d->R_tot[k]=R;
    }
}

static void region_shift_g(Dev *d, const double *dV){
    size_t N=NN(d);
    unsigned char *mark=calloc(N,1); if(!mark) return;
    size_t *q=(size_t*)d->wk[0];
    for(int c=0;c<d->n_cons;c++){
        if(d->cons[c].bc!=0 || fabs(dV[c])<1e-15) continue;
        size_t head=0,tail=0; int type=0;
        for(size_t k=0;k<N;k++) if(d->con_mask[k]==c+1 && d->kind[k]==K_SEMI && !mark[k]){
            mark[k]=(unsigned char)(c+1); q[tail++]=k;
            if(!type){ double D=d->Nd[k]-d->Na[k]; type = D>0.?1:(D<0.?-1:0); }
        }
        if(!type) continue;
        while(head<tail){
            size_t k=q[head++];
            for(int e=d->grp[k];e<d->grp[k+1];e++){
                size_t kb=d->gcol[e];
                if(mark[kb] || d->kind[kb]!=K_SEMI) continue;
                if(d->con_mask[kb] && d->con_mask[kb]!=c+1) continue;
                double D=d->Nd[kb]-d->Na[kb]; if(D*type<=0.) continue;
                double maj = type>0 ? d->n[kb] : d->p[kb];
                if(maj<0.5*fabs(D)) continue;
                mark[kb]=(unsigned char)(c+1); q[tail++]=kb;
            }
        }
        for(size_t t=0;t<tail;t++){
            size_t k=q[t]; if(d->con_mask[k]) continue;
            d->phi[k]+=dV[c]; d->phi_n[k]+=dV[c]; d->phi_p[k]+=dV[c];
        }
    }
    free(mark);
}

static double Emax_g(Dev *d){
    size_t N=NN(d); double Em=0.;
    for(size_t k=0;k<N;k++){
        if(d->kind[k]!=K_SEMI) continue;
        double g[3]; g_grad(d,d->phi,k,g);
        double a=g[0]*g[0]+g[1]*g[1]+g[2]*g[2]; if(a>Em) Em=a;
    }
    return sqrt(Em);
}

/* ======================================================= public API */
static void mark_changed(Dev *d){ d->dirty=1; d->have_solution=0; }

void api3_init(void){
    free_dev(&gdev);
    gdev.n_mats=0; gdev.n_cons=0; gdev.Nx=gdev.Ny=gdev.Nz=0;
    gdev.converged=0; gdev.conv_iter=0; gdev.conv_res=0.;
    gdev.max_iter_p=2000; gdev.max_iter_c=1000; gdev.max_gummel=500;
    gdev.tol_p=1e-9; gdev.tol_c=1e-6; gdev.omega_p=gdev.omega_c=1.;
    gdev.pcg_its=gdev.bicg_its=gdev.newton_its=gdev.gummel_its=0;
    gdev.t_pois=gdev.t_cont=gdev.t_mob=gdev.t_cur=gdev.t_aa=gdev.t_mg=0.;
    { const char *e=getenv("SEMISIM_PCGRTOL"); gdev.pcg_rtol = e? atof(e) : 1e-3;
      const char *f=getenv("SEMISIM_NTOLF");   gdev.ntol_fac = f? atof(f) : 0.1;
      const char *g=getenv("SEMISIM_BRTOL");   gdev.bicg_rtol = g? atof(g) : 1e-4; }
    int nt=1;
#ifdef _OPENMP
    nt=omp_get_max_threads();
#endif
    gdev.nthreads=nt;
    gdev.models=7;                /* BGN + tau(N) + velocity saturation */
    const char *v=getenv("SEMISIM_VERBOSE"); gdev.verbose = v ? atoi(v) : 0;
    const char *a=getenv("SEMISIM_MGA"); if(a) MG_ALPHA=atof(a);
    { const char *m=getenv("SEMISIM_AAM"); if(m){ int q=atoi(m); AA_DEPTH = q<0?0:(q>AA_MAX?AA_MAX:q); } }
    { const char *m=getenv("SEMISIM_CATOL"); if(m) CONT_ATOL=atof(m); }
    { const char *m=getenv("SEMISIM_QFM"); if(m) QF_MARGIN=atof(m); }
}

void api3_set_threads(int n){
#ifdef _OPENMP
    if(n<1) n=omp_get_num_procs();
    omp_set_num_threads(n); gdev.nthreads=n;
#else
    (void)n; gdev.nthreads=1;
#endif
}
void api3_set_verbose(int v){ gdev.verbose=v; }
/* model switches: bit0 band-gap narrowing (Klaassen), bit1 SRH tau(N),
   bit2 field-dependent mobility (velocity saturation), bit3 quantum
   confinement (MLDA), bit4 Lombardi surface mobility, bit5 gate tunnelling */
void api3_set_models(int flags){ gdev.models=flags; gdev.have_solution=0; }

int api3_add_material(int id, double T){
    if(gdev.n_mats>=NM_MAX) return -1;
    mat_init(&gdev.mats[gdev.n_mats],id,T);
    mark_changed(&gdev);
    return gdev.n_mats++;
}

int api3_set_grid_coords(int Nx,int Ny,int Nz,const double *xs,const double *ys,const double *zs){
    for(int i=1;i<Nx;i++) if(!(xs[i]>xs[i-1])) return 0;
    for(int j=1;j<Ny;j++) if(!(ys[j]>ys[j-1])) return 0;
    for(int l=1;l<Nz;l++) if(!(zs[l]>zs[l-1])) return 0;
    if(!alloc_dev(&gdev,Nx,Ny,Nz)) return 0;
    for(int i=0;i<Nx;i++) gdev.xs[i]=xs[i];
    for(int j=0;j<Ny;j++) gdev.ys[j]=ys[j];
    for(int l=0;l<Nz;l++) gdev.zs[l]=zs[l];
    for(int i=0;i<Nx;i++) gdev.Lx[i]=fv_len(gdev.xs,i,Nx);
    for(int j=0;j<Ny;j++) gdev.Ly[j]=fv_len(gdev.ys,j,Ny);
    for(int l=0;l<Nz;l++) gdev.Lz[l]=fv_len(gdev.zs,l,Nz);
    return 1;
}

void api3_set_grid(int Nx,int Ny,int Nz,double Lx,double Ly,double Lz){
    if(Nx<2)Nx=2; if(Ny<2)Ny=2; if(Nz<2)Nz=2;
    double *xs=malloc(Nx*sizeof(double)),*ys=malloc(Ny*sizeof(double)),*zs=malloc(Nz*sizeof(double));
    if(xs&&ys&&zs){
        for(int i=0;i<Nx;i++) xs[i]=i*Lx/(Nx-1.);
        for(int j=0;j<Ny;j++) ys[j]=j*Ly/(Ny-1.);
        for(int l=0;l<Nz;l++) zs[l]=l*Lz/(Nz-1.);
        (void)api3_set_grid_coords(Nx,Ny,Nz,xs,ys,zs);
    }
    free(xs); free(ys); free(zs);
}

/* mp: PCG max its (Poisson), mc: BiCGSTAB max its (continuity),
   mg: Gummel its per continuation step, tp: Gummel tolerance (V),
   tc: relative terminal-current tolerance. op/oc: unused (ABI). */
void api3_set_solver(int mp,int mc,int mg,double tp,double tc,double op,double oc){
    gdev.max_iter_p=mp; gdev.max_iter_c=mc; gdev.max_gummel=mg;
    gdev.tol_p=tp; gdev.tol_c=tc; gdev.omega_p=op; gdev.omega_c=oc;
}

static void clip_box(int *x0,int *y0,int *z0,int *x1,int *y1,int *z1){
    if(*x0<0)*x0=0; if(*y0<0)*y0=0; if(*z0<0)*z0=0;
    if(*x1>=gdev.Nx)*x1=gdev.Nx-1; if(*y1>=gdev.Ny)*y1=gdev.Ny-1; if(*z1>=gdev.Nz)*z1=gdev.Nz-1;
}

void api3_set_region(int x0,int y0,int z0,int x1,int y1,int z1,int midx,double Nd,double Na){
    if(!gdev.allocated||midx<0||midx>=gdev.n_mats) return;
    clip_box(&x0,&y0,&z0,&x1,&y1,&z1);
    int Nx=gdev.Nx,Ny=gdev.Ny;
    for(int l=z0;l<=z1;l++) for(int j=y0;j<=y1;j++) for(int i=x0;i<=x1;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        gdev.mat_idx[k]=midx; gdev.Nd[k]=Nd; gdev.Na[k]=Na;
        if(gdev.mats[midx].insulator){ gdev.Nd[k]=gdev.Na[k]=0.; }
    }
    mark_changed(&gdev);
}

/* Volume fixed charge (m^-3, signed, in units of q) added over a box */
void api3_add_fixed_charge(int x0,int y0,int z0,int x1,int y1,int z1,double Nf){
    if(!gdev.allocated) return;
    clip_box(&x0,&y0,&z0,&x1,&y1,&z1);
    for(int l=z0;l<=z1;l++) for(int j=y0;j<=y1;j++) for(int i=x0;i<=x1;i++)
        gdev.Nf[IDX(i,j,l,gdev.Nx,gdev.Ny)]+=Nf;
    mark_changed(&gdev);
}

/* Sheet charge sigma (elementary charges per m^2, signed) on the plane
   axis=0: x=xs[idx] (a->j, b->l); axis=1: y=ys[idx] (a->i, b->l);
   axis=2: z=zs[idx] (a->i, b->j). Stored as Nf = sigma / L_axis[idx]. */
void api3_add_sheet_charge(int axis,int idx,int a0,int a1,int b0,int b1,double sigma){
    if(!gdev.allocated) return;
    int x0,x1,y0,y1,z0,z1; double L;
    if(axis==0){ if(idx<0||idx>=gdev.Nx) return; x0=x1=idx; y0=a0;y1=a1; z0=b0;z1=b1; L=gdev.Lx[idx]; }
    else if(axis==1){ if(idx<0||idx>=gdev.Ny) return; y0=y1=idx; x0=a0;x1=a1; z0=b0;z1=b1; L=gdev.Ly[idx]; }
    else { if(idx<0||idx>=gdev.Nz) return; z0=z1=idx; x0=a0;x1=a1; y0=b0;y1=b1; L=gdev.Lz[idx]; }
    api3_add_fixed_charge(x0,y0,z0,x1,y1,z1,sigma/L);
}

/* Contacts. face 0..5 = domain faces (i/j = in-plane index ranges, v6
   convention); face 6 = volume box (i: x, j: y, k: z ranges).
   bc 0 ohmic, 1 Schottky (phi_m), 2 gate (phi_m, oxide tox (EOT, m),
   fixed oxide charge Nox (q/m^2)). A gate on semiconductor nodes uses the
   oxide Robin condition; a gate on insulator nodes or a box gate is a
   metal electrode (Dirichlet). Returns the contact id. */
int api3_add_contact_ex(int face,int i0,int i1,int j0,int j1,int k0,int k1,
                        int bc,double V,double pm,double tox,double Nox){
    if(gdev.n_cons>=NC_MAX) return -1;
    Con *c=&gdev.cons[gdev.n_cons];
    memset(c,0,sizeof(*c));
    c->face=face<0?0:(face>6?6:face);
    c->i0=i0;c->i1=i1;c->j0=j0;c->j1=j1;c->k0=k0;c->k1=k1;
    c->bc=bc<0?0:(bc>2?2:bc); c->V=V; c->phi_m=pm;
    c->tox = tox>0.? tox : 10e-9;
    c->qox = Nox*Q;
    { int id=gdev.n_cons; free(gdev.cnodes[id]); free(gdev.careas[id]);
      gdev.cnodes[id]=NULL; gdev.careas[id]=NULL; gdev.cnn[id]=0; }
    mark_changed(&gdev);
    return gdev.n_cons++;
}
void api3_add_contact(int face,int i0,int i1,int j0,int j1,int bc,double V,double pm){
    (void)api3_add_contact_ex(face,i0,i1,j0,j1,0,0,bc,V,pm,10e-9,0.);
}
void api3_clear_contacts(void){ gdev.n_cons=0; mark_changed(&gdev); }
/* Layer stack of a face (Robin) gate for tunnelling, semiconductor side
   first: material indices (api3_add_material) and thicknesses (m).
   nl = 0 removes it (no tunnelling through that gate). */
int api3_set_gate_stack(int cid,int nl,const int *mat,const double *t){
    if(cid<0||cid>=gdev.n_cons||nl<0||nl>TUN_MAXL) return 0;
    Con *c=&gdev.cons[cid];
    for(int i=0;i<nl;i++){ if(mat[i]<0||mat[i]>=gdev.n_mats||!(t[i]>0.)) return 0; }
    c->snl=nl;
    for(int i=0;i<nl;i++){ c->smat[i]=mat[i]; c->st[i]=t[i]; }
    mark_changed(&gdev);
    return 1;
}
/* ------------------------------------------------ unstructured meshes
   A device on a general box-method mesh (e.g. triangular prisms): N nodes,
   E directed edges as a symmetric CSR graph (rp[N+1], cj[E]); per directed
   edge its length h (m) and dual face area ar (m^2); per node the control
   volume vol (m^3) and position pos (3N, m).  Replaces the tensor grid
   (api3_set_grid_coords); materials must be added first.  Returns 1 on
   success (0: bad graph, e.g. not symmetric, or out of memory). */
int api3_set_mesh_graph(int N,int E,const int *rp,const int *cj,const double *h,
                        const double *ar,const double *vol,const double *pos){
    int ok=alloc_dev_graph(&gdev,N,E,rp,cj,h,ar,vol,pos);
    /* The aggregation multigrid of a general mesh reduces smooth error
       components less than the tensor multigrid does for the same residual:
       continuity solves stopped at 1e-4 leave enough error to stall the
       Gummel map, 1e-6 converges it in fewer iterations than the tensor mesh
       needs (SEMISIM_BRTOL still overrides). */
    if(ok && !getenv("SEMISIM_BRTOL")) gdev.bicg_rtol=1e-6;
    return ok;
}
int api3_is_graph(void){ return gdev.unstr; }
/* per-node material index, doping and fixed charge (m^-3) of a graph device */
void api3_set_node_data(const int *mat,const double *Nd,const double *Na,const double *Nf){
    Dev *d=&gdev; if(!d->allocated) return;
    size_t N=NN(d);
    for(size_t k=0;k<N;k++){
        int mi=mat[k]; if(mi<0||mi>=d->n_mats) mi=0;
        d->mat_idx[k]=mi;
        int ins=d->mats[mi].insulator;
        d->Nd[k]= ins ? 0. : Nd[k]; d->Na[k]= ins ? 0. : Na[k]; d->Nf[k]= Nf ? Nf[k] : 0.;
    }
    mark_changed(d);
}
/* nodes of contact cid on a graph device; areas = boundary face area of each
   node for a face contact (m^2), ignored for a box contact */
int api3_set_con_nodes(int cid,int n,const int *nodes,const double *areas){
    Dev *d=&gdev;
    if(cid<0||cid>=d->n_cons||n<0) return 0;
    free(d->cnodes[cid]); free(d->careas[cid]); d->cnodes[cid]=NULL; d->careas[cid]=NULL; d->cnn[cid]=0;
    if(n==0) return 1;
    d->cnodes[cid]=malloc((size_t)n*sizeof(int)); d->careas[cid]=malloc((size_t)n*sizeof(double));
    if(!d->cnodes[cid]||!d->careas[cid]) return 0;
    for(int q=0;q<n;q++){
        if(nodes[q]<0||(size_t)nodes[q]>=NN(d)) return 0;
        d->cnodes[cid][q]=nodes[q]; d->careas[cid][q]= areas ? areas[q] : 0.;
    }
    d->cnn[cid]=n; mark_changed(d);
    return 1;
}
/* interface walls of a graph device: 6 per node - two walls, each the node's
   distance to it and the range a..b of the node's cell measured from it (m);
   distance >= 0.5 m = no wall (MLDA factor 1, no surface mobility) */
void api3_set_walls(const double *w){
    Dev *d=&gdev; if(!d->allocated||!d->unstr) return;
    memcpy(d->gwall,w,NN(d)*6*sizeof(double)); mark_changed(d);
}

/* number of tunnelling paths of the last prepared device (all or contact cid) */
int api3_get_ntp(int cid){
    if(cid<0) return gdev.n_tp;
    return (cid<gdev.n_cons)? gdev.cons[cid].ntp : 0;
}
/* Test hook: tunnelling fluxes (1/m^2 s; e s->g, e g->s, h s->g, h g->s)
   through a stack (layers from the semiconductor surface to the metal)
   from a semiconductor (material index smat) whose surface is at potential
   psi with quasi-Fermi potentials phin/phip, classical carrier densities at
   the surface (local supply), to a metal at Vg, work function pm (eV), in
   the gauge of the current material table.  Returns the net current
   density into the semiconductor through the gate (A/m^2). */
double api3_tun_eval(int smat,int nl,const int *mat,const double *t,double psi,
                     double phin,double phip,double Vg,double pm,double *out){
    Dev *d=&gdev; double F[4]={0,0,0,0};
    if(smat<0||smat>=d->n_mats||nl<1||nl>TUN_MAXL) return 0.;
    for(int i=0;i<nl;i++) if(mat[i]<0||mat[i]>=d->n_mats) return 0.;
    mat_prepare(d);
    const Mat *ms=&d->mats[smat]; TPath P;
    tun_stack(d,&P,nl,mat,t,psi,Vg,pm);
    double Ec=-P.psis-ms->chi, Ev=Ec-ms->Eg;
    tun_flux(ms,&P,Ec,Ev,-phin-d->Cref,-phip-d->Cref,-Vg-d->Cref,F);
    if(out) for(int i=0;i<4;i++) out[i]=F[i];
    return Q*((F[0]-F[1])-(F[2]-F[3]));
}
/* MLDA density factor of material mi at the cell range [a,b] (m) from a wall */
double api3_mlda_factor(int mi,double a,double b,int hole){
    if(mi<0||mi>=gdev.n_mats) return 1.;
    return mlda_factor(&gdev.mats[mi],a,b,hole);
}
/* Lombardi edge mobility (m^2/Vs) of material mi for bulk mobility mub,
   normal field F (V/m), total doping Nt (m^-3), distance factor Dw */
double api3_lombardi(int mi,double mub,double F,double Nt,double Dw,int hole){
    if(mi<0||mi>=gdev.n_mats) return mub;
    return lombardi(&gdev.mats[mi],mub,F,Nt,Dw,hole);
}
void api3_set_contact_V(int cid,double V){ if(cid>=0&&cid<gdev.n_cons) gdev.cons[cid].V=V; }
int    api3_get_ncon(void){ return gdev.n_cons; }
double api3_get_contact_I(int cid){ return (cid>=0&&cid<gdev.n_cons)?gdev.cons[cid].I:0.; }
double api3_get_contact_A(int cid){ return (cid>=0&&cid<gdev.n_cons)?gdev.cons[cid].A:0.; }
double api3_get_contact_In(int cid){ return (cid>=0&&cid<gdev.n_cons)?gdev.cons[cid].In:0.; }
/* current resolution floor of this terminal (A): |I| below it is round-off */
double api3_get_contact_Ifloor(int cid){ return (cid>=0&&cid<gdev.n_cons)?gdev.cons[cid].Iflr:0.; }
/* direct flux out of the contact (diagnostics: sum over terminals = Kirchhoff error) */
double api3_get_contact_Iraw(int cid){ return (cid>=0&&cid<gdev.n_cons)?gdev.cons[cid].Iraw:0.; }
int    api3_get_aa_depth(void){ return gdev.aa_m < AA_DEPTH ? gdev.aa_m : AA_DEPTH; }
double api3_get_contact_Vsolved(int cid){ return (cid>=0&&cid<gdev.n_cons)?gdev.cons[cid].Vprev:0.; }

static int do_solve(Dev *d, int keep){
    if(!prepare(d)) return 0;
    if(!keep || !d->have_solution) equilibrium(d);
    int ok=run_continuation(d);
    upd_J(d); mob_display(d);
    d->converged=ok; g_prog[4]=0.;
    return ok;
}
/* Cooperative stop (any thread): the running solve returns at the next
   Gummel iteration with the last converged bias point loaded (return 0,
   api3_get_contact_Vsolved gives the voltages reached).  The flag stays
   set until api3_clear_stop(), so every later solve returns at once. */
void api3_request_stop(void){ g_stop=1; }
void api3_clear_stop(void){ g_stop=0; }
int  api3_stop_requested(void){ return g_stop; }
/* out[0..7]: see g_prog */
void api3_get_progress(double *out){ for(int i=0;i<8;i++) out[i]=g_prog[i]; }
/* Full solve from equilibrium */
int api3_solve(void){ return do_solve(&gdev,0); }
/* Continue from the previous solution to the current contact voltages */
int api3_resolve(void){ return do_solve(&gdev,1); }

/* I-V sweep of contact cid. Returns terminal current density I/A (A/m^2,
   positive = current flowing INTO the device through that terminal). */
void api3_iv_sweep(int cid,double Vs,double Ve,int Npts,double *Varr,double *Jarr,int *Carr){
    Dev *d=&gdev;
    if(cid<0||cid>=d->n_cons) return;
    if(Npts<1)Npts=1; if(Npts>400)Npts=400;
    double dV=(Npts>1)?(Ve-Vs)/(Npts-1.):0.;
    double V0=d->cons[cid].V;
    if(!prepare(d)) return;
    if(!d->have_solution) equilibrium(d);
    for(int kk=0;kk<Npts;kk++){
        double V=Vs+kk*dV;
        if(g_stop){ Varr[kk]=NAN; Jarr[kk]=NAN; Carr[kk]=-1; continue; }
        d->cons[cid].V=V;
        int ok=run_continuation(d);
        double A=d->cons[cid].A;
        Varr[kk]=d->cons[cid].Vprev;
        Jarr[kk]= A>0. ? d->cons[cid].I/A : d->cons[cid].I;
        Carr[kk]=ok;
    }
    upd_J(d); mob_display(d); g_prog[4]=0.;
    d->cons[cid].V=V0;
}

/* Sweep with all terminal currents: Iarr/Farr are Npts x n_cons (row-major),
   current (A, into device) and its round-off floor (A) of every contact. */
void api3_iv_sweep2(int cid,double Vs,double Ve,int Npts,double *Varr,
                    double *Iarr,double *Farr,int *Carr){
    Dev *d=&gdev;
    if(cid<0||cid>=d->n_cons) return;
    if(Npts<1)Npts=1; if(Npts>400)Npts=400;
    double dV=(Npts>1)?(Ve-Vs)/(Npts-1.):0.;
    double V0=d->cons[cid].V; int nc=d->n_cons;
    if(!prepare(d)) return;
    if(!d->have_solution) equilibrium(d);
    for(int kk=0;kk<Npts;kk++){
        if(g_stop){ Varr[kk]=NAN; Carr[kk]=-1;
            for(int c=0;c<nc;c++){ Iarr[kk*nc+c]=NAN; Farr[kk*nc+c]=NAN; } continue; }
        d->cons[cid].V=Vs+kk*dV;
        int ok=run_continuation(d);
        Varr[kk]=d->cons[cid].Vprev;
        for(int c=0;c<nc;c++){ Iarr[kk*nc+c]=d->cons[c].I; Farr[kk*nc+c]=d->cons[c].Iflr; }
        Carr[kk]=ok;
    }
    upd_J(d); mob_display(d); g_prog[4]=0.;
    d->cons[cid].V=V0;
}

/* Max |E| (V/m) in semiconductor nodes */
double api3_get_Emax(void){
    Dev *d=&gdev; if(!d->allocated) return 0.;
    if(d->unstr) return Emax_g(d);
    int Nx=d->Nx,Ny=d->Ny,Nz=d->Nz; double Em=0.;
    for(int l=0;l<Nz;l++) for(int j=0;j<Ny;j++) for(int i=0;i<Nx;i++){
        size_t k=IDX(i,j,l,Nx,Ny);
        if(d->kind[k]!=K_SEMI) continue;
        double g2=0.;
        for(int ax=0;ax<3;ax++){
            size_t km,kp; double hm,hp,ar;
            int a=nbr(d,i,j,l,2*ax,&km,&hm,&ar), b=nbr(d,i,j,l,2*ax+1,&kp,&hp,&ar);
            double g=0.;
            if(a&&b) g=(d->phi[kp]-d->phi[km])/(hm+hp);
            else if(b) g=(d->phi[kp]-d->phi[k])/hp;
            else if(a) g=(d->phi[k]-d->phi[km])/hm;
            g2+=g*g;
        }
        if(g2>Em) Em=g2;
    }
    return sqrt(Em);
}

int    api3_get_Nx(void){ return gdev.Nx; }
int    api3_get_Ny(void){ return gdev.Ny; }
int    api3_get_Nz(void){ return gdev.Nz; }
int    api3_get_conv(void){ return gdev.converged; }
int    api3_get_iter(void){ return gdev.conv_iter; }
double api3_get_res(void) { return gdev.conv_res; }
double api3_get_xs(int i) { return gdev.xs[i]; }
double api3_get_ys(int j) { return gdev.ys[j]; }
double api3_get_zs(int l) { return gdev.zs[l]; }
double api3_get_Cref(void){ return gdev.Cref; }
long   api3_get_stat(int w){
    switch(w){ case 0: return gdev.pcg_its; case 1: return gdev.bicg_its;
               case 2: return gdev.newton_its; case 3: return gdev.gummel_its; }
    return 0;
}
/* profile: 0 Poisson, 1 continuity, 2 mobility, 3 currents, 4 Anderson, 5 MG setup (s) */
double api3_get_time(int w){
    switch(w){ case 0: return gdev.t_pois; case 1: return gdev.t_cont; case 2: return gdev.t_mob;
               case 3: return gdev.t_cur; case 4: return gdev.t_aa; case 5: return gdev.t_mg; }
    return 0.;
}
double api3_mat_ni(int mi){ return (mi>=0&&mi<gdev.n_mats)?gdev.mats[mi].ni:0.; }
double api3_mat_Eg(int mi){ return (mi>=0&&mi<gdev.n_mats)?gdev.mats[mi].Eg:0.; }

/* Carrier current density component (A/m^2) at node k: mean of the two
   adjacent face fluxes along axis ax; hole=0 electrons, 1 holes */
static double node_Jc(const Dev *d, size_t k, int ax, int hole){
    if(d->kind[k]!=K_SEMI) return 0.;
    if(d->unstr){ double J[3]; g_node_J(d,k,1+hole,J); return J[ax]; }
    int i,j,l; unidx(d,k,&i,&j,&l);
    double s=0.; int cnt=0; size_t kb; double h,ar,Jn,Jp,Sn,Sp;
    if(nbr(d,i,j,l,2*ax,&kb,&h,&ar) && d->kind[kb]==K_SEMI){
        face_flux(d,kb,k,h,ar,&Jn,&Jp,&Sn,&Sp); s+=(hole?Jp:Jn)/ar; cnt++; }
    if(nbr(d,i,j,l,2*ax+1,&kb,&h,&ar) && d->kind[kb]==K_SEMI){
        face_flux(d,k,kb,h,ar,&Jn,&Jp,&Sn,&Sp); s+=(hole?Jp:Jn)/ar; cnt++; }
    return cnt? s/cnt : 0.;
}

/* Derived quantity at node k (display gauge: vacuum level = -phi) */
static double node_val(const Dev *d, int which, size_t k){
    const Mat *m=&d->mats[d->mat_idx[k]]; int kd=d->kind[k];
    switch(which){
    case 0: return d->phi[k];
    case 1: return d->n[k];
    case 2: return d->p[k];
    case 3: return kd==K_METAL? NAN : -d->phi[k]-m->chi-d->bgn[k]*m->VT;
    case 4: return kd==K_METAL? NAN : -d->phi[k]-m->chi-m->Eg+d->bgn[k]*m->VT;
    case 5: return kd==K_INS? NAN : -d->phi_n[k]-d->Cref;
    case 6: return kd==K_INS? NAN : -d->phi_p[k]-d->Cref;
    case 7: return d->Jx[k];
    case 8: return d->Jy[k];
    case 9: return d->Jz[k];
    case 10: return d->mu_n[k];
    case 11: return d->mu_p[k];
    case 12: return d->R_tot[k];
    case 13: return (double)kd;
    case 14: return d->Nd[k]-d->Na[k];
    case 15: return d->phi_n[k];
    case 16: return d->phi_p[k];
    case 17: return d->Nf[k];
    case 18: return (double)d->con_mask[k];
    case 19: return d->bgn[k]*m->VT*2.;           /* band-gap narrowing (eV) */
    case 20: return node_Jc(d,k,0,0);             /* Jn_x */
    case 21: return node_Jc(d,k,1,0);             /* Jn_y */
    case 22: return node_Jc(d,k,2,0);             /* Jn_z */
    case 23: return node_Jc(d,k,0,1);             /* Jp_x */
    case 24: return node_Jc(d,k,1,1);             /* Jp_y */
    case 25: return node_Jc(d,k,2,1);             /* Jp_z */
    case 26: return d->kind[k]==K_SEMI ? niK(d,k) : 0.;  /* effective ni */
    case 27: return d->kind[k]==K_SEMI ? d->tg[k] : 0.;  /* gate tunnelling J (A/m^2) */
    case 28: return d->kind[k]==K_SEMI ? -d->qcn[k]*m->VT : NAN;  /* quantum potential, electrons (eV) */
    case 29: return d->kind[k]==K_SEMI ? -d->qcp[k]*m->VT : NAN;  /* quantum potential, holes (eV) */
    case 30: return d->kind[k]==K_SEMI ? d->qd[k] : NAN;          /* distance to interface (m) */
    }
    return 0.;
}
/* which = 0..30: phi,n,p,Ec,Ev,Efn,Efp,Jx,Jy,Jz,mun,mup,R,kind,Nnet,phin,phip,Nf,
   contact, dEg_bgn, Jn_x,Jn_y,Jn_z, Jp_x,Jp_y,Jp_z, ni_eff, J_tunnel,
   Lambda_n, Lambda_p (quantum potentials, eV), interface distance (SI units) */
void api3_bulk_copy(int which,double *out){
    Dev *d=&gdev; if(!d->allocated) return;
    size_t N=NN(d);
    #pragma omp parallel for schedule(static) if(N>PAR_MIN)
    for(size_t k=0;k<N;k++) out[k]=node_val(d,which,k);
}
double api3_get_phi(int k){ return node_val(&gdev,0,(size_t)k); }
double api3_get_n(int k)  { return node_val(&gdev,1,(size_t)k); }
double api3_get_p(int k)  { return node_val(&gdev,2,(size_t)k); }
double api3_get_Ec(int k) { return node_val(&gdev,3,(size_t)k); }
double api3_get_Ev(int k) { return node_val(&gdev,4,(size_t)k); }
double api3_get_Efn(int k){ return node_val(&gdev,5,(size_t)k); }
double api3_get_Efp(int k){ return node_val(&gdev,6,(size_t)k); }
double api3_get_Jx(int k) { return node_val(&gdev,7,(size_t)k); }
double api3_get_Jy(int k) { return node_val(&gdev,8,(size_t)k); }
double api3_get_Jz(int k) { return node_val(&gdev,9,(size_t)k); }
double api3_get_mun(int k){ return node_val(&gdev,10,(size_t)k); }
double api3_get_mup(int k){ return node_val(&gdev,11,(size_t)k); }
double api3_get_R(int k)  { return node_val(&gdev,12,(size_t)k); }
