#!/usr/bin/env python3
"""SemiSim 3D — 3-D drift-diffusion semiconductor device simulator.
Rectangular (graded tensor) or triangular-prism meshes (semimesh.py);
dense nodes at junctions, interfaces and gates, coarse elsewhere."""

import sys,os,ctypes,threading,time,traceback
# Use all available CPU cores for OpenMP-parallel C solver, unless user
# already set OMP_NUM_THREADS in their environment.
if "OMP_NUM_THREADS" not in os.environ:
    try:
        n_cpu = len(os.sched_getaffinity(0))   # respects taskset/cpusets
    except (AttributeError, OSError):
        n_cpu = os.cpu_count() or 1
    os.environ["OMP_NUM_THREADS"] = str(max(1, n_cpu))
# Short OpenMP spin-wait before sleeping: ~3% slower on an idle machine, but
# 2-3x faster when other programs compete for the cores (long spinning
# threads then wait on descheduled ones at every barrier).
os.environ.setdefault("GOMP_SPINCOUNT", "30000")
import tkinter as tk
from tkinter import ttk,messagebox,filedialog
import numpy as np
import matplotlib; matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D   # noqa
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
try:
    import semimesh as SM                  # triangular-prism mesher (same folder)
except ImportError as e:
    print(f"semimesh.py must be in the same folder as main.py ({e})"); sys.exit(1)

# ══════════════════════════════════════════════════════════════
#  C LIBRARY
# ══════════════════════════════════════════════════════════════
_SO = os.environ.get("SEMISIM_SO") or os.path.join(os.path.dirname(os.path.abspath(__file__)),"semiconductor_core.so")
try: _lib=ctypes.CDLL(_SO)
except OSError as e: print(f"Cannot load .so: {e}"); sys.exit(1)

def _f(fn,res,*args):
    f=getattr(_lib,fn); f.restype=res; f.argtypes=list(args); return f
D=ctypes.c_double; I=ctypes.c_int; V=None; DP=ctypes.POINTER(D)

api_init    = _f("api3_init",V)
api_addmat  = _f("api3_add_material",I,I,D)
api_grid    = _f("api3_set_grid",V,I,I,I,D,D,D)          # uniform fallback
api_gridc   = _f("api3_set_grid_coords",I,I,I,I,DP,DP,DP) # → 1=success, 0=alloc fail
api_solver  = _f("api3_set_solver",V,I,I,I,D,D,D,D)
api_region  = _f("api3_set_region",V,I,I,I,I,I,I,I,D,D)
api_addcon  = _f("api3_add_contact",V,I,I,I,I,I,I,D,D)
# face,i0,i1,j0,j1,k0,k1,bc,V,phi_m,tox(m),Nox(q/m^2) -> contact id
api_addconx = _f("api3_add_contact_ex",I,I,I,I,I,I,I,I,I,D,D,D,D)
api_clrcon  = _f("api3_clear_contacts",V)
api_sheet   = _f("api3_add_sheet_charge",V,I,I,I,I,I,I,D)  # axis,idx,a0,a1,b0,b1,sigma(q/m^2)
api_models  = _f("api3_set_models",V,I)                    # bit0 BGN, bit1 tau(N), bit2 v-sat
api_solve   = _f("api3_solve",I)
api_resolve = _f("api3_resolve",I)
api_setV    = _f("api3_set_contact_V",V,I,D)
api_conI    = _f("api3_get_contact_I",D,I)      # terminal current (A, into device)
api_conA    = _f("api3_get_contact_A",D,I)      # contact area (m^2)
api_conIn   = _f("api3_get_contact_In",D,I)     # electron component (A)
api_conF    = _f("api3_get_contact_Ifloor",D,I) # round-off resolution floor (A)
api_stat    = _f("api3_get_stat",ctypes.c_long,I)
api_Nx=_f("api3_get_Nx",I); api_Ny=_f("api3_get_Ny",I); api_Nz=_f("api3_get_Nz",I)
api_conv=_f("api3_get_conv",I); api_iter=_f("api3_get_iter",I); api_res=_f("api3_get_res",D)
api_xs=_f("api3_get_xs",D,I); api_ys=_f("api3_get_ys",D,I); api_zs=_f("api3_get_zs",D,I)
_Gd={nm:_f(f"api3_get_{nm}",D,I) for nm in
     ["phi","n","p","Ec","Ev","Efn","Efp","Jx","Jy","Jz","mun","mup","R"]}
api_mat_ni=_f("api3_mat_ni",D,I); api_mat_Eg=_f("api3_mat_Eg",D,I)
api_bulk=_f("api3_bulk_copy",V,I,DP)  # which, out
IP=ctypes.POINTER(I)
api_Emax=_f("api3_get_Emax",D)        # max E-field (V/m) for breakdown
_IV_N=400
_V_arr=(D*_IV_N)(); _J_arr=(D*_IV_N)(); _C_arr=(I*_IV_N)()
api_iv=_f("api3_iv_sweep",V,I,D,D,I,DP,DP,ctypes.POINTER(I))
# sweep returning every terminal current + its resolution floor (Npts x ncon)
api_iv2=_f("api3_iv_sweep2",V,I,D,D,I,DP,DP,DP,ctypes.POINTER(I))
api_init()
BYTES_PER_NODE=760   # v7 core incl. multigrid hierarchy and depth-12 Anderson history (depth drops if RAM is short)

# ══════════════════════════════════════════════════════════════
#  NON-UNIFORM GRID GENERATOR
# ══════════════════════════════════════════════════════════════
def detect_junctions(regs, Lx, Ly, Lz, doping_ratio_thresh=10.0):
    """Auto-detect junction planes (positions where doping or material
    changes abruptly) from a list of Reg objects.

    Returns (dense_x, dense_y, dense_z) — sorted lists of nm positions
    along each axis where mesh refinement is physically warranted.

    Logic:
      For each pair of regions (A, B) that overlap or touch in 2D
      (i.e. share a face plane along one axis) AND have either
      (a) different material, (b) different doping type (n vs p),
      or (c) doping ratio > thresh — record the shared boundary
      coordinate as a junction point.

    This excludes single-region domain centers (like y=100 in a
    bulk-uniform diode) from the dense list, fixing the spurious
    middle-of-bulk refinement seen in the old fallback.
    """
    if not regs: return [], [], []

    def signed_doping(r):
        """Net doping with sign: + for n, − for p, 0 for intrinsic."""
        if r.doping <= 1e11: return 0.0
        return r.doping if r.dtype == "n" else -r.doping

    def regions_overlap_2d(rA, rB, axis):
        """Do regions A and B overlap on the two axes OTHER than `axis`?"""
        if axis == "x":
            return (rA.y0 < rB.y1 and rB.y0 < rA.y1 and
                    rA.z0 < rB.z1 and rB.z0 < rA.z1)
        if axis == "y":
            return (rA.x0 < rB.x1 and rB.x0 < rA.x1 and
                    rA.z0 < rB.z1 and rB.z0 < rA.z1)
        return (rA.x0 < rB.x1 and rB.x0 < rA.x1 and
                rA.y0 < rB.y1 and rB.y0 < rA.y1)

    def physically_different(rA, rB):
        """Decide if a junction at the shared plane is physically meaningful."""
        if rA.mat != rB.mat:
            return True   # heterojunction (e.g. AlGaAs/GaAs)
        dA = signed_doping(rA); dB = signed_doping(rB)
        if (dA > 0 and dB < 0) or (dA < 0 and dB > 0):
            return True   # pn junction
        # high-low junction (e.g. n+/n-): same sign but ratio >= thresh
        if dA != 0 and dB != 0:
            r = max(abs(dA), abs(dB)) / min(abs(dA), abs(dB))
            if r >= doping_ratio_thresh:
                return True
        # intrinsic-to-doped boundary
        if (dA == 0) != (dB == 0):
            return True
        return False

    junx, juny, junz = set(), set(), set()
    n = len(regs)
    for i in range(n):
        for j in range(i+1, n):
            A, B = regs[i], regs[j]
            if not physically_different(A, B):
                continue
            # For each axis, find boundary planes where A and B meet AND
            # they overlap on the OTHER two axes. We need to check edges of
            # BOTH regions because an inner small region's edges (e.g. a
            # source island inside a body) are only on that small region's
            # x0/x1 — not on the surrounding region.
            for x_val in [A.x0, A.x1, B.x0, B.x1]:
                if 0 < x_val < Lx and regions_overlap_2d(A, B, "x"):
                    # x_val is a junction if it lies inside the *other*
                    # region's x extent (so the two regions touch/overlap
                    # there along x).
                    other_lo, other_hi = (B.x0, B.x1) if x_val in (A.x0, A.x1) else (A.x0, A.x1)
                    if other_lo <= x_val <= other_hi:
                        junx.add(round(x_val, 1))
            for y_val in [A.y0, A.y1, B.y0, B.y1]:
                if 0 < y_val < Ly and regions_overlap_2d(A, B, "y"):
                    other_lo, other_hi = (B.y0, B.y1) if y_val in (A.y0, A.y1) else (A.y0, A.y1)
                    if other_lo <= y_val <= other_hi:
                        juny.add(round(y_val, 1))
            for z_val in [A.z0, A.z1, B.z0, B.z1]:
                if 0 < z_val < Lz and regions_overlap_2d(A, B, "z"):
                    other_lo, other_hi = (B.z0, B.z1) if z_val in (A.z0, A.z1) else (A.z0, A.z1)
                    if other_lo <= z_val <= other_hi:
                        junz.add(round(z_val, 1))

    # Drop the domain-edge points (refinement at boundary just shifts nodes
    # without resolving anything new — boundary already has nodes).
    junx = sorted(p for p in junx if 0 < p < Lx)
    juny = sorted(p for p in juny if 0 < p < Ly)
    junz = sorted(p for p in junz if 0 < p < Lz)
    return junx, juny, junz


def graded_grid(L_nm, N, dense_nm, ratio=8.0, sigma_nm=None):
    """
    Generate N non-uniform nodes in [0, L_nm] nm.
    Nodes are ~ratio× denser near positions in dense_nm than far away.
    Uses CDF inversion of a Gaussian density function.

    sigma_nm: width (nm) of refinement region around each dense point.
              If None, auto-chosen as min(half-distance to nearest dense pt, L/20).
              Smaller sigma → tighter refinement, more uniform spacing in bulk.

    Returns coordinate array in SI metres.
    """
    L = L_nm * 1e-9
    if not dense_nm or N < 3:
        return np.linspace(0, L, N)
    # Auto-pick sigma per dense point if not given.
    # Each peak gets sigma = min(half-distance to nearest other peak, L/20).
    # This guarantees the refinement region stays LOCAL — no leak into middle.
    pts = sorted(set(dense_nm))
    if sigma_nm is None:
        sigmas = []
        for p in pts:
            others = [q for q in pts if q != p]
            if not others:
                d_min = L_nm * 0.10   # standalone peak: sigma = 10% of domain
            else:
                d_min = min(abs(p - q) for q in others) * 0.5   # half-distance
            sigmas.append(min(d_min, L_nm * 0.05) * 1e-9)        # cap at 5% of L
    else:
        sigmas = [sigma_nm * 1e-9] * len(pts)
    # Floor sigma at ~1.5x the uniform spacing so we don't get duplicate nodes
    sigmas = [max(s, L / N * 1.5) for s in sigmas]

    x_test = np.linspace(0, L, 20000)
    density = np.ones(20000)
    for d, sig in zip(pts, sigmas):
        density += (ratio - 1.0) * np.exp(-0.5 * ((x_test - d*1e-9) / sig)**2)
    cdf = np.cumsum(density); cdf -= cdf[0]; cdf /= cdf[-1]
    xs = np.interp(np.linspace(0, 1, N), cdf, x_test)
    xs[0] = 0.; xs[-1] = L
    for i in range(1, N):       # guarantee strict monotone
        if xs[i] <= xs[i-1]: xs[i] = xs[i-1] + 1e-13
    return xs

def uniform_grid(L_nm, N):
    return np.linspace(0., L_nm*1e-9, N)

def span_idx(coords_m, a_m, b_m):
    """Inclusive node-index range (i0,i1) for the physical HALF-OPEN interval
    [a, b) on a (possibly graded) monotonic grid.

    Using the same half-open rule for regions and contacts guarantees that
      * two regions or contacts sharing a boundary never share a node,
      * a contact drawn over a region only ever touches that region's nodes
        (the old inclusive/nearest-node rules let a BJT base contact reach one
        node into the collector - a base-collector short),
      * placement is by physical coordinate, so graded meshes are handled
        correctly (the old contact rule assumed a uniform mesh).
    An interval that reaches the far domain edge includes the last node; an
    interval too thin to contain any node gets the node nearest its centre."""
    c=np.asarray(coords_m); N=len(c); L=c[-1]-c[0]; tol=1e-9*max(L,1e-30)
    if a_m<=c[0]+tol: i0=0
    else: i0=int(np.searchsorted(c, a_m-tol, side="left"))
    if b_m>=c[-1]-tol: i1=N-1
    else: i1=int(np.searchsorted(c, b_m-tol, side="left"))-1
    i0=max(0,min(i0,N-1)); i1=max(-1,min(i1,N-1))
    if i1<i0:
        k=int(np.argmin(np.abs(c-0.5*(a_m+b_m)))); i0=i1=k
    return i0,i1

def coord_to_idx(xs, x_nm):
    """Return index closest to x_nm (nm) in xs (m) array."""
    x = x_nm * 1e-9
    return int(np.argmin(np.abs(xs - x)))

# ══════════════════════════════════════════════════════════════
#  THEME
# ══════════════════════════════════════════════════════════════
MATS  = {"Si":0,"Ge":1,"GaAs":2,"InP":3,"4H-SiC":4,"GaN":5,"AlGaAs":6,"AlGaN":7,
         "SiO2":8,"Al2O3":9,"HfO2":10,"Si3N4":11}
INSULATORS = {"SiO2","Al2O3","HfO2","Si3N4"}
EPS_R = {"SiO2":3.9,"Al2O3":9.0,"HfO2":22.0,"Si3N4":7.5}   # insulators: EOT of gate stacks
MCOLS = {"Si":("#3b9de8","#e85050"),"Ge":("#2ec08a","#f07030"),
         "GaAs":("#a855f7","#e8305a"),"InP":("#06b6d4","#e05090"),
         "4H-SiC":("#22c55e","#eab308"),"GaN":("#38bdf8","#f97316"),
         "AlGaAs":("#8b5cf6","#f43f5e"),"AlGaN":("#67e8f9","#fb923c"),
         "SiO2":("#cbd5e1","#cbd5e1"),"Al2O3":("#e2e8f0","#e2e8f0"),
         "HfO2":("#fcd34d","#fcd34d"),"Si3N4":("#d6d3d1","#d6d3d1")}
C={"bg":"#05090f","panel":"#0b1422","border":"#182436",
   "accent":"#00a8f0","green":"#00cc80","yellow":"#f0c020",
   "red":"#ff3050","purple":"#b060ff","orange":"#ff7820",
   "text":"#dae8fb","sub":"#5a8aaa",
   "Ec":"#38bdf8","Ev":"#f87171","Efn":"#86efac","Efp":"#fbbf24"}
TCOLS=["#00a8f0","#00cc80","#b060ff","#f0c020","#ff7820","#ff3050",
       "#06b6d4","#8b5cf6","#10b981","#ec4899","#14b8a6","#f59e0b","#ef4444","#6366f1"]
plt.rcParams.update({
    "figure.facecolor":C["bg"],"axes.facecolor":C["bg"],"axes.edgecolor":C["border"],
    "axes.labelcolor":C["sub"],"xtick.color":C["sub"],"ytick.color":C["sub"],
    "text.color":C["text"],"grid.color":C["border"],"grid.alpha":0.45,
    "legend.facecolor":C["panel"],"legend.edgecolor":C["border"],
    "legend.labelcolor":C["text"],"legend.fontsize":7})

def sax(ax,xl="",yl="",t="",log=False):
    ax.set_xlabel(xl,fontsize=8); ax.set_ylabel(yl,fontsize=8)
    ax.set_title(t,fontsize=9,fontweight="bold",pad=4); ax.grid(True,lw=0.35,ls="--")
    if log: ax.set_yscale("log")
def sax3(ax,xl="",yl="",zl="",t=""):
    try: ax.locator_params(nbins=4)
    except Exception: pass
    ax.set_xlabel(xl,fontsize=7,color=C["sub"]); ax.set_ylabel(yl,fontsize=7,color=C["sub"])
    ax.set_zlabel(zl,fontsize=7,color=C["sub"]); ax.set_title(t,fontsize=9,fontweight="bold",color=C["text"])
    ax.tick_params(labelsize=6,colors=C["sub"]); ax.set_facecolor(C["bg"])
    ax.xaxis.pane.fill=ax.yaxis.pane.fill=ax.zaxis.pane.fill=False
def mfig(p,fs=(11,5.5)):
    """Figure + canvas + matplotlib toolbar.  The toolbar is packed FIRST at
    the bottom so a short window squeezes the plot, never the toolbar; the
    small initial size lets the layout, not the figure, decide the size."""
    fig=Figure(figsize=(5,3.5),facecolor=C["bg"])
    cv=FigureCanvasTkAgg(fig,p)
    tf=tk.Frame(p,bg="#020609"); tf.pack(side=tk.BOTTOM,fill=tk.X)
    try:
        tb=NavigationToolbar2Tk(cv,tf,pack_toolbar=False); tb.pack(side=tk.LEFT,fill=tk.X)
    except TypeError:
        tb=NavigationToolbar2Tk(cv,tf)
    tb.config(bg="#020609"); tb.update()
    cv.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,expand=True)
    return fig,cv

# ══════════════════════════════════════════════════════════════
#  DATA CLASSES
# ══════════════════════════════════════════════════════════════
class Reg:
    """Region: a box (nm), optionally cut to a cross-section outline (shape,
    semimesh.Shape: the region is its box intersected with the outline, or
    its box minus the outline for a hole shape).  Later regions override."""
    def __init__(self,mat,x0,y0,z0,x1,y1,z1,dtype,doping,label="",shape=None):
        self.mat=mat; self.x0=x0;self.y0=y0;self.z0=z0
        self.x1=x1;self.y1=y1;self.z1=z1; self.dtype=dtype
        self.doping=doping; self.label=label or f"{mat} {dtype} {doping:.0e}"
        self.shape=SM.Shape.from_dict(shape) if shape is not None else None
class Con:
    FNAME=["X-Min","X-Max","Y-Min","Y-Max","Z-Min","Z-Max","Box"]
    BCNAME=["Ohmic","Schottky","Gate"]
    # Axis labels: for each face, what do i% and j% correspond to? (Box: i→X j→Y k→Z)
    FAXES={0:("Y","Z"),1:("Y","Z"),2:("X","Z"),3:("X","Z"),4:("X","Y"),5:("X","Y"),6:("X","Y")}
    def __init__(self,face,i0p,i1p,j0p,j1p,bc,V=0.,pm=4.05,label="",
                 tox=10.,nox=0.,k0p=0.,k1p=1.,net="",stack=None,shape=None):
        self.face=face; self.i0p=i0p;self.i1p=i1p;self.j0p=j0p;self.j1p=j1p
        self.k0p=k0p; self.k1p=k1p
        # Box contacts only: cross-section outline (semimesh.Shape) cutting the
        # box - e.g. a gate electrode = box minus the outline of its gate stack
        self.shape=SM.Shape.from_dict(shape) if (shape is not None and face==6) else None
        self.bc=bc; self.V=V; self.pm=pm
        self.tox=tox      # gate oxide thickness, nm (EOT, SiO2 permittivity)
        self.nox=nox      # fixed oxide charge at the interface, q/cm^2
        # gate stack of a face gate, semiconductor side first: [(material, nm)];
        # None/empty = a single SiO2 layer of thickness tox.  With a stack the
        # gate capacitance uses its EOT and tunnelling sees every layer.
        self.stack=[(str(m),float(t)) for m,t in stack] if stack else None
        self.label=label or f"{Con.FNAME[face]} {Con.BCNAME[bc]} {V:.1f}V"
        # net: contacts sharing a net name are one electrical node (wired
        # together outside the simulated cell) - they always carry the same
        # voltage, are swept together and their currents add up
        self.net=net or ""
    def to_idx(self,Nx,Ny,Nz,xs=None,ys=None,zs=None):
        """Fraction coords (0..1) -> inclusive node ranges (i0,i1,j0,j1).
        Pass the actual node coordinates (xs,ys,zs) whenever a grid exists:
        the contact is then placed by PHYSICAL position with the same half-open
        rule as regions (see span_idx). Without coordinates a uniform mesh is
        assumed - only acceptable for display before meshing."""
        if xs is not None and ys is not None and zs is not None:
            ax={0:(ys,zs),1:(ys,zs),2:(xs,zs),3:(xs,zs),4:(xs,ys),5:(xs,ys),6:(xs,ys)}[self.face]
            c1=np.asarray(ax[0]); c2=np.asarray(ax[1])
            i0,i1=span_idx(c1, c1[0]+self.i0p*(c1[-1]-c1[0]), c1[0]+self.i1p*(c1[-1]-c1[0]))
            j0,j1=span_idx(c2, c2[0]+self.j0p*(c2[-1]-c2[0]), c2[0]+self.j1p*(c2[-1]-c2[0]))
            return i0,i1,j0,j1
        if self.face in(0,1):   fs1=max(Ny-1,1); fs2=max(Nz-1,1)
        elif self.face in(2,3): fs1=max(Nx-1,1); fs2=max(Nz-1,1)
        else:                   fs1=max(Nx-1,1); fs2=max(Ny-1,1)
        i0=int(round(self.i0p*fs1)); i1=int(round(self.i1p*fs1))
        j0=int(round(self.j0p*fs2)); j1=int(round(self.j1p*fs2))
        i0=max(0,min(i0,fs1)); i1=max(i0,min(i1,fs1))
        j0=max(0,min(j0,fs2)); j1=max(j0,min(j1,fs2))
        return i0,i1,j0,j1
    def to_idx6(self,xs,ys,zs):
        """(i0,i1,j0,j1,k0,k1) for the C core; k-range only used by Box."""
        i0,i1,j0,j1=self.to_idx(len(xs),len(ys),len(zs),xs,ys,zs)
        if self.face!=6: return i0,i1,j0,j1,0,0
        c=np.asarray(zs)
        k0,k1=span_idx(c, c[0]+self.k0p*(c[-1]-c[0]), c[0]+self.k1p*(c[-1]-c[0]))
        return i0,i1,j0,j1,k0,k1
    def layers(self):
        """Gate stack actually used: [(material, nm)], semiconductor side first."""
        return list(self.stack) if self.stack else [("SiO2",float(self.tox))]
    def eot(self):
        """Equivalent oxide thickness (nm) of the gate stack."""
        if not self.stack: return float(self.tox)
        return sum(t*3.9/EPS_R.get(m,3.9) for m,t in self.stack)
    @staticmethod
    def stack_str(stack):
        return ", ".join(f"{m} {t:g}" for m,t in stack) if stack else ""
    @staticmethod
    def parse_stack(txt):
        """'SiO2 0.5, HfO2 2.8' (nm, semiconductor side first) -> [(mat, nm)] or None."""
        txt=(txt or "").strip()
        if not txt: return None
        out=[]
        for part in txt.replace(";",",").split(","):
            w=part.replace(":"," ").replace("="," ").split()
            if not w: continue
            if len(w)!=2 or w[0] not in INSULATORS: raise ValueError(f"bad layer '{part.strip()}'")
            t=float(w[1])
            if not t>0: raise ValueError(f"thickness must be > 0 in '{part.strip()}'")
            out.append((w[0],t))
        if len(out)>6: raise ValueError("at most 6 layers")
        return out or None
    def axes_hint(self):
        """Return string showing what i%/j% axes mean for this face."""
        a=Con.FAXES[self.face]
        return f"i→{a[0]} j→{a[1]}" + (" k→Z" if self.face==6 else "")

def con_full(face,bc,V,lbl="",pm=4.05,tox=10.,nox=0.,net="",stack=None):
    return Con(face,0.,1.,0.,1.,bc,V,pm,lbl,tox,nox,net=net,stack=stack)
def con_part(face,i0p,i1p,j0p,j1p,bc,V,lbl="",pm=4.05,tox=10.,nox=0.,net="",stack=None):
    return Con(face,i0p,i1p,j0p,j1p,bc,V,pm,lbl,tox,nox,net=net,stack=stack)
def con_box(x0,y0,z0,x1,y1,z1,L,bc,V,lbl="",pm=4.05,net="",shape=None):
    """Buried electrode occupying the box (nm) in a domain of size L=(Lx,Ly,Lz).
    A Gate box is a metal/poly electrode (e.g. inside an oxide-lined trench);
    an Ohmic box on semiconductor is a buried contact.  shape: optional
    cross-section outline cutting the box (see Reg)."""
    Lx,Ly,Lz=L
    return Con(6,x0/Lx,x1/Lx,y0/Ly,y1/Ly,bc,V,pm,lbl,10.,0.,z0/Lz,z1/Lz,net=net,shape=shape)

def obj_dict(o):
    """vars() of a Reg/Con/Sheet, JSON-ready (shapes as dicts)."""
    d=dict(vars(o))
    if d.get("shape") is not None: d["shape"]=d["shape"].to_dict()
    return d

class Sheet:
    """Fixed interface sheet charge sigma (q/cm^2, signed) on the plane
    <axis> = pos (nm), limited to the ranges a0..a1, b0..b1 (nm) of the two
    other axes taken in x,y,z order. Represents e.g. the AlGaN/GaN
    polarisation charge that induces the 2DEG."""
    def __init__(self,axis,pos,a0,a1,b0,b1,sigma,label=""):
        self.axis=axis; self.pos=pos; self.a0=a0; self.a1=a1; self.b0=b0; self.b1=b1
        self.sigma=sigma; self.label=label or f"sheet {sigma:+.2e} q/cm2 @{axis}={pos:.0f}nm"

# ══════════════════════════════════════════════════════════════
#  DEVICE TEMPLATES — with dense-point specs
# ══════════════════════════════════════════════════════════════
# Each template returns a dict with:
#   Lx,Ly,Lz, Nx,Ny,Nz, regs, cons, ivc, ivS, ivE, ivN, desc
#   dense_x, dense_y, dense_z: lists of nm positions for grid refinement

def T_pn():
    Lx,Ly,Lz=400.,200.,200.; Nx,Ny,Nz=36,20,6
    # junction at x=200nm, contacts at x=0 and x=400
    regs=[Reg("Si",0,0,0,200,200,200,"p",1e17,"p-Si anode"),
          Reg("Si",200,0,0,400,200,200,"n",1e17,"n-Si cathode")]
    cons=[con_full(0,0,0.,"Anode"),con_full(1,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-0.8,ivE=0.8,ivN=33,
                dense_x=[200.],dense_y=[],dense_z=[],
                desc="Si p-n diode — junction dense at x=200nm")

def T_pin():
    Lx,Ly,Lz=600.,200.,200.; Nx,Ny,Nz=42,20,6
    regs=[Reg("Si",0,0,0,100,200,200,"p",2e18,"p+ emitter"),
          Reg("Si",100,0,0,500,200,200,"n",1e14,"i-region"),
          Reg("Si",500,0,0,600,200,200,"n",2e18,"n+ collector")]
    cons=[con_full(0,0,0.,"Anode"),con_full(1,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-1.,ivE=1.,ivN=41,
                dense_x=[100.,500.],dense_y=[],dense_z=[],
                desc="Si PIN diode — dense at both p+/i and i/n+ junctions")

def T_schottky():
    Lx,Ly,Lz=500.,200.,200.; Nx,Ny,Nz=38,20,6
    regs=[Reg("Si",0,0,0,500,200,200,"n",1e16,"n-Si"),
          Reg("Si",450,0,0,500,200,200,"n",1e18,"n+ ohmic back")]
    cons=[con_full(0,1,0.,"Schottky Metal",4.72),con_full(1,0,0.,"Ohmic Back")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-1.5,ivE=0.7,ivN=45,
                dense_x=[0.,450.],dense_y=[],dense_z=[],
                desc="Ti/n-Si Schottky diode — dense near metal interface")

def _bjt(pol):
    """Vertical bipolar transistor cross-section (y = depth, top surface y=0):
       n+ emitter (top-left) | p+ extrinsic base under the base contact
       (top-right) | p base | n collector | n+ buried layer, collector contact
       at the bottom. The base contact sits on a p+ region beside the emitter,
       as in a planar process, so injected electrons cross the INTRINSIC base
       vertically instead of recombining at the base metal."""
    Lx,Ly,Lz=2000.,1600.,200.; Nx,Ny,Nz=52,70,4
    e,b = ("n","p") if pol=="npn" else ("p","n")
    s_=1. if pol=="npn" else -1.
    regs=[Reg("Si",0,0,0,Lx,Ly,Lz,        e,1e16,f"{e} collector"),
          Reg("Si",0,1300,0,Lx,Ly,Lz,     e,1e19,f"{e}+ buried layer"),
          Reg("Si",0,0,0,Lx,450,Lz,       b,1e17,f"{b} base (0.25um)"),
          Reg("Si",0,0,0,900,200,Lz,      e,2e19,f"{e}+ emitter"),
          Reg("Si",1400,0,0,Lx,200,Lz,    b,1e19,f"{b}+ extrinsic base")]
    cons=[con_part(2,0.,900/Lx,0.,1.,0,0.0,"Emitter"),
          con_part(2,1400/Lx,1.,0.,1.,0,0.7*s_,"Base"),
          con_full(3,0,2.0*s_,"Collector")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=0.35*s_,ivE=0.80*s_,ivN=19,ivout="Collector",
                dense_x=[900.,1400.],dense_y=[200.,450.,1300.],dense_z=[],
                desc=f"{pol.upper()} BJT - vertical, p+ extrinsic base; sweep V_BE (Gummel plot of all terminals)")
def T_npn(): return _bjt("npn")
def T_pnp(): return _bjt("pnp")

def T_nmos():
    """Planar n-MOSFET (y = depth, top surface y = 0 like every template).
       n+ source/drain 100 nm deep at the surface, n+ poly gate (phi_m 4.1 eV)
       over 7 nm SiO2 between them (L = 320 nm), p-body 3e17 with the bulk
       contact at the bottom.  V_th ~ 0.6 V (0.57 V classical + ~25 mV from
       quantum confinement).  Surface (Lombardi) mobility lowers the on-current
       by ~30 % against bulk mobility; gate tunnelling through 7 nm is
       negligible.  Default V_G = V_D = 1.5 V; the sweep is the output curve
       I_D(V_D) (set sweep = Gate for I_D(V_G))."""
    Lx,Ly,Lz=500.,300.,200.; Nx,Ny,Nz=42,30,6
    regs=[Reg("Si",0,0,0,500,300,200,"p",3e17,"p-body bulk"),
          Reg("Si",0,0,0,90,100,200,"n",1e20,"n+ source"),     # top-left
          Reg("Si",410,0,0,500,100,200,"n",1e20,"n+ drain")]   # top-right
    # Source, gate and drain on the top face (face 2 = Y-Min):
    cons=[con_part(2,0.,0.18,0.,1.,0,0.,"Source"),       # over n+ source
          con_part(2,0.82,1.,0.,1.,0,1.5,"Drain"),       # over n+ drain
          con_part(2,0.18,0.82,0.,1.,2,1.5,"Gate",4.1,tox=7.),  # n+ poly, 7 nm oxide
          con_full(3,0,0.,"Bulk")]                        # bottom
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=0.,ivE=3.,ivN=31,wl=200./320.,models=63,
                dense_x=[90.,410.],dense_y=[0.,100.],dense_z=[],
                desc="N-MOSFET - 7nm oxide gate, Vth~0.5V; source/gate/drain on top, bulk at bottom")

def T_pmos():
    """Planar p-MOSFET: complementary to the NMOS template (p+ poly gate,
    phi_m 5.2 eV, n-well 3e17, V_th ~ -0.52 V).  Holes lose about half of
    their bulk mobility to the surface (Lombardi), so the PMOS delivers
    ~45 % of the NMOS current at the same |V|."""
    Lx,Ly,Lz=500.,300.,200.; Nx,Ny,Nz=42,30,6
    regs=[Reg("Si",0,0,0,500,300,200,"n",3e17,"n-well bulk"),
          Reg("Si",0,0,0,90,100,200,"p",1e20,"p+ source"),
          Reg("Si",410,0,0,500,100,200,"p",1e20,"p+ drain")]
    cons=[con_part(2,0.,0.18,0.,1.,0,0.,"Source"),
          con_part(2,0.82,1.,0.,1.,0,-1.5,"Drain"),
          con_part(2,0.18,0.82,0.,1.,2,-1.5,"Gate",5.2,tox=7.),  # p+ poly
          con_full(3,0,0.,"N-well")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=-3.,ivE=0.,ivN=31,wl=200./320.,models=63,
                dense_x=[90.,410.],dense_y=[0.,100.],dense_z=[],
                desc="P-MOSFET - 7nm oxide, p+ poly gate; source/gate/drain on top, n-well at bottom")

def T_vdmos():
    """Vertical Double-diffused MOSFET (VDMOS) — power switch
    Truly VERTICAL device with proper JFET window between p-body cells.

    Block stacking — KEY: p-body has a GAP in the middle so n- drift
    extends up to the surface there (the "JFET window"). Current path:
    n+ source → lateral channel under gate → JFET window → drift → drain.

      Layer 1: n- drift everywhere (background)
      Layer 2,3: p-body LEFT and RIGHT cells (with JFET gap in middle)
      Layer 4: n+ substrate at bottom (drain)
      Layer 5,6: n+ source islands at top (carve into p-body)
      Layer 7,8: p+ body tap (carve into p-body, next to source)

    Cross-section (X horizontal, Y is depth, top at y=0):

       y=0  ┌[n+src]─[p+tap]─┤ JFET ├─[p+tap]─[n+src]┐  ← top
            │  p-body     ←gate→     ←gate→  p-body  │
            │             │ DRIFT │              │   │
            │~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~│
            │       n- drift bulk (long for high BV) │
            ├────────────────────────────────────────┤
       y=Ly │      n+ substrate                      │  ← bottom (drain)
            └────────────────────────────────────────┘
    """
    Lx,Ly,Lz=4000.,12000.,400.; Nx,Ny,Nz=44,72,4
    body_dep = 1000.   # p-body cell depth
    src_dep  = 200.    # n+ source / p+ tap depth
    sub_top  = 11000.
    src_w    = 800.    # n+ source island width (from outer edge inward)
    tap_w    = 300.    # p+ body tap width (next to source)
    body_w   = src_w + tap_w + 200.   # p-body cell extends slightly past tap
    # JFET window centered: x in [Lx/2 - jw/2, Lx/2 + jw/2]
    jfet_hw  = 600.    # JFET window half-width

    bx0_L = 0.;        bx1_L = Lx/2. - jfet_hw   # left p-body block
    bx0_R = Lx/2. + jfet_hw; bx1_R = Lx          # right p-body block

    regs=[
        # Layer 1: background — n- drift
        Reg("Si",0,0,0,Lx,Ly,Lz, "n",5e15, "n- drift"),
        # Layer 2: p-body LEFT (with JFET gap on right side)
        Reg("Si",bx0_L,0,0,bx1_L,body_dep,Lz, "p",2e17, "p-body L"),
        # Layer 3: p-body RIGHT
        Reg("Si",bx0_R,0,0,bx1_R,body_dep,Lz, "p",2e17, "p-body R"),
        # Layer 4: n+ substrate at bottom
        Reg("Si",0,sub_top,0,Lx,Ly,Lz, "n",5e19, "n+ substrate"),
        # Layer 5,6: p+ body tap at the outer edge of each cell (shorted to source)
        Reg("Si",0,0,0,tap_w,src_dep,Lz, "p",1e19, "p+ tap L"),
        Reg("Si",Lx-tap_w,0,0,Lx,src_dep,Lz, "p",1e19, "p+ tap R"),
        # Layer 7,8: n+ source next to the channel (tap | source | channel)
        Reg("Si",tap_w,0,0,tap_w+src_w,src_dep,Lz, "n",1e20, "n+ source L"),
        Reg("Si",Lx-tap_w-src_w,0,0,Lx-tap_w,src_dep,Lz, "n",1e20, "n+ source R"),
    ]
    # Source contact spans p+ tap + n+ source on each side
    src_frac = (src_w + tap_w) / Lx
    # Gate: overlaps the source edge by 200 nm, covers the channel (source edge
    # to p-body edge) and extends 200 nm over the JFET window (accumulation)
    gate_L_x0 = tap_w + src_w - 200.
    gate_L_x1 = Lx/2. - jfet_hw + 200.
    gate_R_x0 = Lx/2. + jfet_hw - 200.
    gate_R_x1 = Lx - tap_w - src_w + 200.
    cons=[
        con_part(2, 0., src_frac,                       0.,1., 0, 0., "Source-L"),
        con_part(2, 1.-src_frac, 1.,                    0.,1., 0, 0., "Source-R"),
        con_part(2, gate_L_x0/Lx, gate_L_x1/Lx,         0.,1., 2, 10., "Gate-L", 4.1, tox=50.),
        con_part(2, gate_R_x0/Lx, gate_R_x1/Lx,         0.,1., 2, 10., "Gate-R", 4.1, tox=50.),
        con_full(3, 0, 5., "Drain"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=4,ivS=0.,ivE=20.,ivN=21,
                dense_x=[tap_w, tap_w+src_w, bx1_L, bx0_R, Lx-tap_w-src_w, Lx-tap_w],
                dense_y=[0., src_dep, body_dep, sub_top],
                dense_z=[],
                desc="VDMOS - vertical; tap|source|channel, 50nm gate oxide, JFET window")

def T_ldmos():
    """Lateral Double-diffused MOSFET (LDMOS) — RF & audio power amplifier
    LATERAL device: source/gate/drain all on top, but with extended n- drift
    region between gate and drain to support high V_DS.

    Block stacking:
      Layer 1: p- substrate (background)
      Layer 2: n- drift extended drain region (CARVES top half)
      Layer 3: p-body (CARVES at source side, must overlap drift edge)
      Layer 4: n+ source island
      Layer 5: n+ drain
      Layer 6: p+ body tap (under source, shorted to source metal)

    Cross-section:
       y=0   ┌──[src][bodytap]─[gate]──[ n- drift ]──[drain]┐ ← top
             │  n+    p+      |gateox|                  n+   │
             │       p-body                                  │ ← p-body wraps source
             │       p- substrate (long bulk)                │
             └───────────────────────────────────────────────┘ ← bottom (substrate tap)
    """
    Lx,Ly,Lz=2400.,800.,300.; Nx,Ny,Nz=44,30,4
    src_w   = 250.   # source island width
    tap_w0,tap_w1 = 250., 500.  # body tap x range
    body_w  = 800.   # p-body extent in x (covers source/tap/channel)
    drift_x = 700.   # n-drift starts here (must overlap p-body for proper junction)
    drain_w = 200.   # n+ drain width
    surf_d  = 200.   # depth of surface S/D regions
    drift_d = 250.   # drift depth
    body_d  = 350.   # p-body depth

    regs=[
        # Layer 1: p- substrate (background)
        Reg("Si",0,0,0,Lx,Ly,Lz, "p",1e16, "p- substrate"),
        # Layer 2: n- drift in top region from drift_x to drain
        Reg("Si",drift_x,0,0,Lx,drift_d,Lz, "n",2e16, "n- drift (ext drain)"),
        # Layer 3: p-body (CARVES drift in source-side overlap region)
        Reg("Si",0,0,0,body_w,body_d,Lz, "p",5e17, "p-body"),
        # Layer 4: p+ body tap at the outer edge
        Reg("Si",0,0,0,src_w,surf_d,Lz, "p",1e19, "p+ body tap"),
        # Layer 5: n+ source next to the channel (tap | source | channel)
        Reg("Si",tap_w0,0,0,tap_w1,surf_d,Lz, "n",5e19, "n+ source"),
        # Layer 6: n+ drain at far right
        Reg("Si",Lx-drain_w,0,0,Lx,surf_d,Lz, "n",5e19, "n+ drain"),
    ]
    cons=[
        # Source metal shorts n+ source and p+ body tap (one wide contact)
        con_part(2, 0., tap_w1/Lx,                  0.,1.,0,0.,    "Source"),
        # Gate: overlaps the source edge, covers the channel, extends over the drift
        con_part(2, (tap_w1-100.)/Lx, (body_w+100.)/Lx, 0.,1.,2,5.,  "Gate", 4.1, tox=25.),
        # Drain: top, far right
        con_part(2, (Lx-drain_w)/Lx, 1.,            0.,1.,0,15.,   "Drain"),
        # Substrate body contact at bottom (often grounded for RF LDMOS)
        con_full(3, 0, 0., "Substrate"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=2,ivS=0.,ivE=20.,ivN=21,bv_thr=0.05,
                dense_x=[src_w,tap_w1,body_w,Lx-drain_w],
                dense_y=[surf_d,body_d,drift_d],
                dense_z=[],
                desc="LDMOS — lateral power FET with extended n- drift, ~20V class")

def T_igbt():
    """Planar field-stop IGBT, 600 V class (one 6 um cell: two half gates).
    MOS gate over a p-body/n- drift/n field-stop/p+ collector stack: channel
    electrons flow down the n- drift and forward-bias the p+ collector
    junction, which injects holes -> the drift is flooded with an e-h plasma
    (conductivity modulation), which is why an IGBT beats a MOSFET of the same
    rating in on-state voltage.
      n- drift 1e14 cm^-3 x 30 um + 2 um n field stop (1e16): trapezoidal
      field, BV_ideal ~ E_c W - qNW^2/2eps ~ 680 V.
      JFET window 2 um with an n JFET implant (1e16): without it the lightly
      doped drift between the p-bodies is pinched off by their built-in
      depletion and the channel current has nowhere to go.
    Default bias is the on-state (V_GE = 15 V, V_CE = 2 V); the sweep is the
    on-state output curve with the ~0.7 V collector-junction knee.  For the
    blocking state set Gate = 0 V and sweep the Collector up to ~600 V."""
    Lx,Ly,Lz=6000.,34000.,400.; Nx,Ny,Nz=48,90,4
    body_dep, src_dep, jfet_dep = 1500., 250., 2000.
    fs_top, coll_top = 31000., 33000.
    src_w, tap_w, jfet_hw = 800., 300., 1000.
    bx1_L = Lx/2. - jfet_hw; bx0_R = Lx/2. + jfet_hw
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz, "n",1e14, "n- drift (600V)"),
        Reg("Si",bx1_L-500,0,0,bx0_R+500,jfet_dep,Lz, "n",1e16, "n JFET implant"),
        Reg("Si",0,0,0,bx1_L,body_dep,Lz, "p",2e17, "p-body L"),
        Reg("Si",bx0_R,0,0,Lx,body_dep,Lz, "p",2e17, "p-body R"),
        Reg("Si",0,fs_top,0,Lx,coll_top,Lz, "n",1e16, "n field stop"),
        Reg("Si",0,coll_top,0,Lx,Ly,Lz, "p",1e18, "p+ collector"),
        Reg("Si",0,0,0,tap_w,src_dep,Lz, "p",1e19, "p+ tap L"),
        Reg("Si",Lx-tap_w,0,0,Lx,src_dep,Lz, "p",1e19, "p+ tap R"),
        Reg("Si",tap_w,0,0,tap_w+src_w,src_dep,Lz, "n",1e20, "n+ emitter L"),
        Reg("Si",Lx-tap_w-src_w,0,0,Lx-tap_w,src_dep,Lz, "n",1e20, "n+ emitter R"),
    ]
    sf=(src_w+tap_w)/Lx
    cons=[
        con_part(2, 0., sf,                                  0.,1., 0, 0., "Emitter-L"),
        con_part(2, 1.-sf, 1.,                               0.,1., 0, 0., "Emitter-R"),
        con_part(2, (tap_w+src_w-100)/Lx, (bx1_L+200)/Lx,     0.,1., 2, 15., "Gate-L", 4.1, tox=100.),
        con_part(2, (bx0_R-200)/Lx, (Lx-tap_w-src_w+100)/Lx, 0.,1., 2, 15., "Gate-R", 4.1, tox=100.),
        con_full(3, 0, 2., "Collector"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=4,ivS=0.,ivE=5.,ivN=26,
                dense_x=[tap_w, tap_w+src_w, bx1_L, bx0_R, Lx-tap_w-src_w, Lx-tap_w],
                dense_y=[0., src_dep, body_dep, jfet_dep, fs_top, coll_top],
                dense_z=[],
                desc="Planar field-stop IGBT 600V - JFET implant, conductivity modulation; on-state sweep")

def T_solar():
    """Si Solar Cell — vertical PV with light entering from top
    n+/p/p+ stack with front grid contact and full-area back contact.

    Block stacking:
      Layer 1: p-base absorber (background)
      Layer 2: n+ emitter (top thin layer)
      Layer 3: p+ BSF (back-surface field, bottom thin layer)
    """
    Lx,Ly,Lz=300.,3000.,300.; Nx,Ny,Nz=18,42,4
    em_dep   = 80.    # n+ emitter depth (~80nm typical)
    bsf_top  = 2900.  # p+ BSF at the back

    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz, "p",1e16, "p-base absorber"),
        Reg("Si",0,0,0,Lx,em_dep,Lz, "n",5e18, "n+ front emitter"),
        Reg("Si",0,bsf_top,0,Lx,Ly,Lz, "p",1e18, "p+ BSF"),
    ]
    cons=[
        # Front contact: small finger on top (real cells have grid contacts)
        con_part(2, 0.30, 0.70,  0.,1., 0, 0., "Front (n+)"),
        # Back contact: full back face (Al BSF contact)
        con_full(3, 0, 0., "Back (p+)"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=-0.1,ivE=0.7,ivN=33,
                dense_x=[Lx/2.],dense_y=[0.,em_dep,bsf_top,Ly],dense_z=[],
                desc="Si solar cell — vertical, light enters top, n+/p/p+ stack")

def T_sic_pn():
    """4H-SiC PiN power diode — high voltage rectifier
    Layered p+/i (n-)/n+ structure for high BV."""
    Lx,Ly,Lz=400.,8000.,400.; Nx,Ny,Nz=20,42,4
    pp_dep   = 500.    # p+ anode top region
    drift_to = 7000.   # n- drift bottom edge
    regs=[
        Reg("4H-SiC",0,0,0,Lx,Ly,Lz, "n",5e15, "n- drift"),
        Reg("4H-SiC",0,0,0,Lx,pp_dep,Lz, "p",1e18, "p+ anode"),
        Reg("4H-SiC",0,drift_to,0,Lx,Ly,Lz, "n",1e19, "n+ substrate"),
    ]
    cons=[
        con_full(2, 0, 0., "Anode"),
        con_full(3, 0, 0., "Cathode"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-1200.,ivE=4.,ivN=41,bv_thr=0.05,
                dense_x=[Lx/2.],dense_y=[0.,pp_dep,drift_to,Ly],dense_z=[],
                desc="4H-SiC PiN power diode — vertical, ~1.2kV class")

def T_gan_hemt():
    """AlGaN/GaN HEMT (depletion mode) - high-frequency RF power FET.
    The 2DEG is induced by the POLARISATION sheet charge at the AlGaN/GaN
    interface (+1.2e13 q/cm^2 for Al0.25Ga0.75N; surface donors are assumed to
    compensate the top-surface polarisation charge), with an UNDOPED 25 nm
    AlGaN barrier - no artificial channel doping, so the 2DEG keeps the
    undoped-GaN mobility.  Ni/Au Schottky gate (phi_m 5.1 eV -> phi_B ~1.3 eV)."""
    Lx,Ly,Lz=1200.,400.,200.; Nx,Ny,Nz=44,44,4
    bar=25.; s_w=200.; d_x=Lx-200.; cdep=60.
    regs=[
        Reg("GaN",0,0,0,Lx,Ly,Lz, "p",1e15, "GaN buffer (SI)"),
        Reg("AlGaN",0,0,0,Lx,bar,Lz, "n",1e10, "AlGaN barrier (undoped)"),
        Reg("GaN",0,0,0,s_w,cdep,Lz, "n",5e19, "n+ source ohmic"),
        Reg("GaN",d_x,0,0,Lx,cdep,Lz, "n",5e19, "n+ drain ohmic"),
    ]
    cons=[
        con_part(2, 0., s_w/Lx,                0.,1., 0, 0., "Source"),
        con_part(2, d_x/Lx, 1.,                0.,1., 0, 5., "Drain"),
        con_part(2, 0.35, 0.50,                0.,1., 1, 0., "Gate (Ni)", 5.1),
    ]
    sheets=[Sheet("y",bar,0.,Lx,0.,Lz,1.2e13,"AlGaN/GaN polarisation")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,sheets=sheets,
                ivc=1,ivS=0.,ivE=10.,ivN=21,
                dense_x=[s_w,d_x,0.35*Lx,0.50*Lx],
                dense_y=[0.,bar,cdep],dense_z=[],
                desc="AlGaN/GaN HEMT - polarisation 2DEG (1.2e13 q/cm2), undoped barrier, Schottky gate")
def T_trench_sic_mos():
    """4H-SiC trench MOSFET, 1.2 kV class.
    SiO2-lined trench (60 nm) with a buried poly gate; vertical channel along
    the sidewall; p+ shield under the trench bottom protects the gate oxide
    from the drift-region field (modern SiC trench feature).
    The shield sits only under the trench bottom, and an n current-spreading
    layer (8e16) keeps the gap between p-body and shield open - without it
    the body/shield built-in depletion (~3 V in SiC) pinches the channel exit.
    Default: on-state (V_GS = 15 V, V_DS = 2 V); sweep = output curve.
    Note: the SiC/SiO2 channel mobility (~20-40 cm^2/Vs in real devices) is
    NOT degraded here - bulk mobility is used, so on-current is optimistic."""
    Lx,Ly,Lz=4000.,12000.,400.; Nx,Ny,Nz=60,90,4
    tw=400.; tx0,tx1=Lx/2-tw,Lx/2+tw; td=1500.; tox=60.
    body_d=1200.; src_d=250.; sub_top=11000.; src_w=800.; tap_w=600.
    regs=[
        Reg("4H-SiC",0,0,0,Lx,Ly,Lz, "n",8e15, "n- drift"),
        Reg("4H-SiC",0,sub_top,0,Lx,Ly,Lz, "n",1e19, "n+ substrate"),
        Reg("4H-SiC",0,body_d,0,Lx,td+700,Lz, "n",8e16, "n current-spreading layer"),
        Reg("4H-SiC",0,0,0,Lx,body_d,Lz, "p",2e17, "p-body"),
        Reg("4H-SiC",tx0,td,0,tx1,td+300,Lz, "p",5e18, "p+ trench-bottom shield"),
        Reg("4H-SiC",tx0-src_w,0,0,tx0,src_d,Lz, "n",1e20, "n+ source L"),
        Reg("4H-SiC",tx1,0,0,tx1+src_w,src_d,Lz, "n",1e20, "n+ source R"),
        Reg("4H-SiC",0,0,0,tap_w,src_d,Lz, "p",5e19, "p+ body tap L"),
        Reg("4H-SiC",Lx-tap_w,0,0,Lx,src_d,Lz, "p",5e19, "p+ body tap R"),
        Reg("SiO2",tx0,0,0,tx1,td,Lz, "i",0., "trench oxide"),
    ]
    L=(Lx,Ly,Lz)
    cons=[
        con_part(2, 0., tx0/Lx,           0.,1., 0, 0., "Source-L"),
        con_part(2, tx1/Lx, 1.,           0.,1., 0, 0., "Source-R"),
        con_box(tx0+tox,0.,0.,tx1-tox,td-tox,Lz, L, 2, 15., "Trench gate (n+ poly)", 4.1),
        con_full(3, 0, 2., "Drain"),
        # the shield is tied to the source at the stripe ends in a real die;
        # in this cross-section that tie is a buried contact inside it
        con_box(tx0+100.,td+80.,0.,tx1-100.,td+220.,Lz, L, 0, 0., "Shield (to source)"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=3,ivS=0.,ivE=5.,ivN=21,
                dense_x=[tap_w,tx0-src_w,tx0,tx0+tox,tx1-tox,tx1,tx1+src_w,Lx-tap_w],
                dense_y=[0.,src_d,body_d,td-tox,td,td+300,td+700,sub_top],
                dense_z=[],
                desc="4H-SiC trench MOSFET 1.2kV - oxide-lined trench, buried gate, p+ shield")
def T_trench_si_mos():
    """Si trench MOSFET, 30 V class (DC-DC converter switch).
    The trench is SiO2-lined (50 nm) and filled with an n+ poly gate electrode
    (buried Gate contact); the channel forms VERTICALLY along the trench
    sidewall in the p-body, between the n+ source and the n- drift."""
    Lx,Ly,Lz=2400.,3500.,300.; Nx,Ny,Nz=60,64,4
    tw=250.; tx0,tx1=Lx/2-tw,Lx/2+tw; td=900.; tox=50.
    body_d=700.; src_d=200.; sub_top=3000.; src_w=500.; tap_w=350.
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz, "n",5e15, "n- epi"),
        Reg("Si",0,sub_top,0,Lx,Ly,Lz, "n",5e19, "n+ substrate"),
        Reg("Si",0,0,0,Lx,body_d,Lz, "p",1.5e17, "p-body"),
        Reg("Si",tx0-src_w,0,0,tx0,src_d,Lz, "n",5e19, "n+ source L"),
        Reg("Si",tx1,0,0,tx1+src_w,src_d,Lz, "n",5e19, "n+ source R"),
        Reg("Si",0,0,0,tap_w,src_d,Lz, "p",1e19, "p+ body tap L"),
        Reg("Si",Lx-tap_w,0,0,Lx,src_d,Lz, "p",1e19, "p+ body tap R"),
        Reg("SiO2",tx0,0,0,tx1,td,Lz, "i",0., "trench oxide"),
    ]
    L=(Lx,Ly,Lz)
    cons=[
        con_part(2, 0., tx0/Lx,           0.,1., 0, 0., "Source-L"),
        con_part(2, tx1/Lx, 1.,           0.,1., 0, 0., "Source-R"),
        con_box(tx0+tox,0.,0.,tx1-tox,td-tox,Lz, L, 2, 10., "Trench gate (n+ poly)", 4.1),
        con_full(3, 0, 1., "Drain"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=3,ivS=0.,ivE=3.,ivN=16,
                dense_x=[tap_w,tx0-src_w,tx0,tx0+tox,tx1-tox,tx1,tx1+src_w,Lx-tap_w],
                dense_y=[0.,src_d,body_d,td-tox,td,sub_top],
                dense_z=[],
                desc="Si trench MOSFET 30V - SiO2-lined trench, buried poly gate, vertical channel")
def T_sic_jbs():
    """4H-SiC Junction Barrier Schottky (JBS) Diode — 1.2kV class
    Schottky cathode top + interspersed p+ stripes that pinch off the
    Schottky region under reverse bias (lower leakage than pure Schottky).

    Block stacking:
      Layer 1: n- drift (background)
      Layer 2: n+ substrate at bottom
      Layer 3: p+ stripes at top (carve into drift to form JBS structure)
    """
    Lx,Ly,Lz=4000.,8000.,500.; Nx,Ny,Nz=36,48,4
    sub_top = 7200.
    pp_d = 400.    # p+ stripe depth
    pp_w = 400.    # p+ stripe width
    pp_x_centers = [600., 2000., 3400.]  # 3 p+ stripes evenly spaced

    regs=[
        Reg("4H-SiC",0,0,0,Lx,Ly,Lz, "n",5e15, "n- drift"),
        Reg("4H-SiC",0,sub_top,0,Lx,Ly,Lz, "n",5e18, "n+ substrate"),
    ]
    # p+ JBS stripes — each carved into drift at top
    for cx in pp_x_centers:
        regs.append(Reg("4H-SiC", cx-pp_w/2, 0, 0, cx+pp_w/2, pp_d, Lz,
                        "p", 1e18, f"p+ JBS stripe @{cx:.0f}nm"))

    cons=[
        # Top: Schottky metal (Ti, phi_m ~4.5 eV) covers entire top
        # The p+ stripes form ohmic-like junctions with Schottky metal
        con_full(2, 1, 0., "Schottky Anode", 4.5),
        con_full(3, 0, 0., "Cathode"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-1200.,ivE=2.,ivN=31,bv_thr=0.05,
                dense_x=[c for c in pp_x_centers] + [c-pp_w/2 for c in pp_x_centers]
                       + [c+pp_w/2 for c in pp_x_centers],
                dense_y=[0.,pp_d,sub_top],
                dense_z=[],
                desc="4H-SiC JBS diode — 1.2kV, Schottky + p+ JBS stripes")

def T_led_gaas():
    Lx,Ly,Lz=500.,200.,200.; Nx,Ny,Nz=40,20,6
    regs=[Reg("GaAs",0,0,0,150,200,200,"p",3e18,"p+-GaAs"),
          Reg("GaAs",150,0,0,350,200,200,"n",2e16,"active i-region"),
          Reg("GaAs",350,0,0,500,200,200,"n",3e18,"n+-GaAs")]
    cons=[con_full(0,0,0.,"Anode"),con_full(1,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-1.,ivE=1.5,ivN=51,
                dense_x=[150.,350.],dense_y=[],dense_z=[],
                desc="GaAs p+/i/n+ LED — Eg=1.42eV  λ≈870nm")

def T_hbt():
    """AlGaAs/GaAs HBT - vertical npn with a wide-bandgap AlGaAs emitter.
       The emitter/base valence-band step blocks hole back-injection, so the
       base can be doped far ABOVE the emitter (4e19 vs 5e17): a thin,
       low-resistance base with high gain - the defining HBT trade.
       Layout (y down): AlGaAs emitter mesa on the left, 80 nm p+ GaAs base,
       n- collector, n+ sub-collector (collector contact at the bottom); the
       base contact sits on a p+ extrinsic-base ledge 200 nm beside the mesa
       (a contact on top of the intrinsic base would swallow the injected
       electrons).  V_BE = 1.35 V (GaAs turn-on ~1.2-1.3 V) -> J_C of order
       1e3 A/cm^2; the sweep is the output curve I_C(V_CE).
       The E-B junction is ABRUPT: its 0.24 eV conduction-band spike throttles
       electron injection (drift-diffusion has no thermionic-emission model),
       so only dEv ~ 0.12 eV of the gap step helps and beta ~ 10.  Production
       HBTs grade the E-B composition over ~30 nm to remove the spike."""
    Lx,Ly,Lz=1200.,1000.,200.; Nx,Ny,Nz=40,60,4
    regs=[Reg("GaAs",0,0,0,Lx,Ly,Lz,"n",5e15,"GaAs collector"),
          Reg("GaAs",0,700,0,Lx,Ly,Lz,"n",1e19,"n+ sub-collector"),
          Reg("GaAs",0,100,0,Lx,180,Lz,"p",4e19,"p+ GaAs base (C-doped)"),
          Reg("AlGaAs",0,0,0,600,100,Lz,"n",5e17,"AlGaAs emitter mesa"),
          Reg("GaAs",600,0,0,Lx,100,Lz,"p",4e19,"p+ extrinsic base")]
    cons=[con_part(2,0.,500/Lx,0.,1.,0,0.0,"Emitter"),
          con_part(2,800/Lx,1.,0.,1.,0,1.35,"Base"),
          con_full(3,0,2.0,"Collector")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=2,ivS=0.,ivE=3.,ivN=31,
                dense_x=[500.,600.,800.],dense_y=[0.,100.,180.,700.],dense_z=[],
                desc="AlGaAs/GaAs HBT - vertical, 5e17 wide-gap emitter over a 4e19 base (abrupt E-B)")

def T_jfet():
    """n-channel JFET, dual gate (top + bottom p+ gates tied together).
    Channel 2e16 cm^-3, 800 nm thick between the gate junctions: pinch-off
    voltage V_p = qN a^2/2eps ~ 2.5 V, so V_GS(off) ~ -(V_p - V_bi) ~ -1.6 V.
    Default V_GS = -0.5 V (channel partly depleted), V_DS = 3 V (saturated:
    the channel pinches at the drain end); sweep = output curve."""
    Lx,Ly,Lz=3000.,1100.,200.; Nx,Ny,Nz=60,44,4
    g0,g1=500.,2500.; gb,gt=150.,950.
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,"n",2e16,"n-channel"),
        Reg("Si",g0,gt,0,g1,Ly,Lz,"p",5e18,"p+ gate top"),
        Reg("Si",g0,0,0,g1,gb,Lz,"p",5e18,"p+ gate bot"),
    ]
    cons=[con_full(0,0,0.,"Source"),
          con_full(1,0,3.,"Drain"),
          con_part(3,g0/Lx,g1/Lx,0.,1.,0,-0.5,"Gate-Top"),
          con_part(2,g0/Lx,g1/Lx,0.,1.,0,-0.5,"Gate-Bot")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=0.,ivE=5.,ivN=26,
                dense_x=[g0,g1],dense_y=[gb,gt],dense_z=[],
                desc="n-JFET - 0.8 um channel between tied p+ gates, V_GS(off) ~ -1.6 V")

# ══════════════════════════════════════════════════════════════
#  REAL-WORLD DEVICE LIBRARY
#
#  These templates are modelled on recognisable commercial parts. Layer
#  thicknesses and dopings are DERIVED, not guessed: for every blocking
#  layer the drift doping and thickness come from the one-sided abrupt
#  junction relations
#         BV = eps*Ec^2 / (2*q*Nd)        Wdrift = 2*BV / Ec
#  with Ec = 3.0e5 V/cm (Si), 2.5e6 V/cm (4H-SiC), 3.3e6 V/cm (GaN),
#  so each structure is self-consistent with its voltage class.
#
#  IMPORTANT / HONEST CAVEATS
#   * Real manufacturers do not publish internal doping profiles. These are
#     physically-consistent representative structures for the device CLASS,
#     not reverse-engineered copies of a specific die.
#   * A simulated unit cell (a few um wide) stands in for a die containing
#     millions of such cells in parallel, so absolute currents are per-cm^2,
#     not per-part amperes.
#   * There is NO impact-ionisation model, so avalanche breakdown voltage is
#     not predicted from first principles; the GUI's breakdown marker is a
#     heuristic knee detector on the I-V curve. Blocking-state fields and
#     depletion widths ARE computed correctly and can be compared with Ec.
#   * Real Zener/tunnel behaviour needs band-to-band tunnelling, also absent.
# ══════════════════════════════════════════════════════════════

def T_rw_diode_1n4148():
    """Si small-signal switching diode, 100 V class (1N4148-like).
    Vertical p+/n/n+ mesa. Drift: Nd=4.3e15 cm^-3, 6 um -> BV_ideal ~ 100 V."""
    Lx,Ly,Lz = 1000., 8000., 1000.; Nx,Ny,Nz = 8, 150, 6
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,        "n",4.3e15,"n drift (100V)"),
        Reg("Si",0,0,0,Lx,400,Lz,       "p",1e19, "p+ anode"),
        Reg("Si",0,6400,0,Lx,Ly,Lz,     "n",1e19, "n+ cathode"),
    ]
    cons=[con_full(2,0,0.,"Anode"), con_full(3,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-120.,ivE=1.0,ivN=41,bv_thr=0.05,
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si switching diode, 100V class (1N4148-like) - p+/n/n+ vertical")

def T_rw_schottky_1n5819():
    """Si Schottky power rectifier, 40 V / 1 A class (1N5819-like).
    Barrier metal phi_m=4.60 eV -> phi_B~0.55 eV (low Vf, higher leakage).
    Epi: Nd=1.2e16 cm^-3, 3 um -> BV_ideal ~ 46 V."""
    Lx,Ly,Lz = 1000., 4500., 1000.; Nx,Ny,Nz = 8, 120, 6
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,        "n",1.2e16,"n epi (40V)"),
        Reg("Si",0,3000,0,Lx,Ly,Lz,     "n",1e19,  "n+ substrate"),
    ]
    cons=[con_full(2,1,0.,"Schottky Anode",4.60), con_full(3,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-50.,ivE=0.6,ivN=41,bv_thr=0.05,
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si Schottky rectifier 40V/1A (1N5819-like) - phi_B~0.55 eV")

def T_rw_zener_5v1():
    """Si Zener voltage reference, 5.1 V class (1N4733A-like).
    Both sides degenerate (p+ 5e19 / n+ 2e19) so the junction is only ~10 nm
    wide and the field reaches >5 MV/cm at 5 V - the tunnelling regime.
    NOTE: no band-to-band tunnelling model, so the 5.1 V knee is NOT
    reproduced. What IS correct and instructive is the ultra-thin, ultra-high
    field depletion layer that makes Zener action possible."""
    Lx,Ly,Lz = 400., 300., 400.; Nx,Ny,Nz = 6, 200, 6
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,      "n",2e19,"n+ cathode side"),
        Reg("Si",0,0,0,Lx,150,Lz,     "p",5e19,"p+ anode side"),
    ]
    cons=[con_full(2,0,0.,"Anode"), con_full(3,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-8.,ivE=0.9,ivN=41,
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si 5.1V Zener (1N4733A-like) - 10nm junction, >5 MV/cm")

def T_rw_pin_photodiode():
    """Si PIN photodiode, visible/NIR (BPW34-like).
    Thick 30 um near-intrinsic absorber gives a wide depletion at low reverse
    bias -> low capacitance, fast response. Run at -5 V (fully depleted)."""
    Lx,Ly,Lz = 2000., 32000., 2000.; Nx,Ny,Nz = 8, 180, 6
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,        "n",5e12, "i (nu) absorber"),
        Reg("Si",0,0,0,Lx,300,Lz,       "p",1e19, "p+ window"),
        Reg("Si",0,30500,0,Lx,Ly,Lz,    "n",1e19, "n+ cathode"),
    ]
    cons=[con_full(2,0,-5.,"Anode (p+)"), con_full(3,0,0.,"Cathode (n+)")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-20.,ivE=0.5,ivN=41,
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si PIN photodiode (BPW34-like) - 30um depleted absorber")

def T_rw_ge_photodiode():
    """Ge PIN photodiode for 1310/1550 nm fibre optics.
    Ge Eg=0.66 eV -> cutoff ~1.88 um, so it covers both telecom windows.
    Thin 5 um absorber because Ge absorbs strongly."""
    Lx,Ly,Lz = 1000., 6200., 1000.; Nx,Ny,Nz = 8, 130, 6
    regs=[
        Reg("Ge",0,0,0,Lx,Ly,Lz,        "n",1e13, "i-Ge absorber"),
        Reg("Ge",0,0,0,Lx,200,Lz,       "p",1e19, "p+ contact"),
        Reg("Ge",0,5200,0,Lx,Ly,Lz,     "n",1e19, "n+ contact"),
    ]
    cons=[con_full(2,0,-2.,"Anode (p+)"), con_full(3,0,0.,"Cathode (n+)")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-10.,ivE=0.4,ivN=41,
                dense_x=[],dense_y=[],dense_z=[],
                desc="Ge PIN photodiode 1310/1550nm telecom - Eg=0.66 eV")

def T_rw_bjt_bc547():
    """Si NPN small-signal BJT, 45 V / beta 200-450 (BC547-like).
    Vertical npn with a laterally-offset p+ base contact, as in a real planar
    process. Base width 0.6 um sets the gain; collector 1e16/4 um sets BV_CEO."""
    Lx,Ly,Lz = 8000., 6400., 1000.; Nx,Ny,Nz = 44, 120, 4
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,          "n",1e16, "n collector"),
        Reg("Si",0,4900,0,Lx,Ly,Lz,       "n",1e19, "n+ subcollector"),
        Reg("Si",0,0,0,Lx,900,Lz,         "p",2e17, "p base (0.6um under emitter)"),
        Reg("Si",0,0,0,3000,300,Lz,       "n",1e20, "n+ emitter"),
        Reg("Si",5000,0,0,Lx,900,Lz,      "p",1e19, "p+ base contact"),
    ]
    cons=[
        con_part(2, 0.,      3000/Lx, 0.,1., 0, 0.0,  "Emitter"),
        con_part(2, 5000/Lx, 1.,      0.,1., 0, 0.7,  "Base"),
        con_full(3, 0, 3.0, "Collector"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=0.,ivE=0.85,ivN=35,ivout="Collector",
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si NPN small-signal BJT 45V (BC547-like) - 0.6um base")

def T_rw_mosfet_irlz44n():
    """Si logic-level power MOSFET, 55 V / 47 A class (IRLZ44N-like).
    One planar VDMOS half-cell with a JFET window. Epi Nd=8e15 cm^-3, 4.5 um
    -> BV_ideal ~ 63 V (margin over the 55 V rating, as in real parts).
    Gate work function 4.1 eV (n+ poly) gives a logic-level Vth."""
    Lx,Ly,Lz = 4000., 8000., 400.; Nx,Ny,Nz = 44, 150, 4
    body_d, src_d, sub_top = 1200., 300., 6000.
    src_w, tap_w, jfet_hw  = 800., 300., 600.
    bx1_L = Lx/2. - jfet_hw; bx0_R = Lx/2. + jfet_hw
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,                    "n",8e15, "n epi (55V)"),
        Reg("Si",0,0,0,bx1_L,body_d,Lz,             "p",2e17, "p-body L"),
        Reg("Si",bx0_R,0,0,Lx,body_d,Lz,            "p",2e17, "p-body R"),
        Reg("Si",0,sub_top,0,Lx,Ly,Lz,              "n",5e19, "n+ substrate"),
        Reg("Si",0,0,0,tap_w,src_d,Lz,              "p",1e19, "p+ tap L"),
        Reg("Si",Lx-tap_w,0,0,Lx,src_d,Lz,          "p",1e19, "p+ tap R"),
        Reg("Si",tap_w,0,0,tap_w+src_w,src_d,Lz,    "n",1e20, "n+ source L"),
        Reg("Si",Lx-tap_w-src_w,0,0,Lx-tap_w,src_d,Lz,"n",1e20,"n+ source R"),
    ]
    sf=(src_w+tap_w)/Lx
    cons=[
        con_part(2, 0., sf,                        0.,1., 0, 0.,  "Source-L"),
        con_part(2, 1.-sf, 1.,                     0.,1., 0, 0.,  "Source-R"),
        con_part(2, (tap_w+src_w-200)/Lx, (bx1_L+200)/Lx, 0.,1., 2, 5., "Gate-L", 4.1, tox=40.),
        con_part(2, (bx0_R-200)/Lx, (Lx-tap_w-src_w+200)/Lx, 0.,1., 2, 5., "Gate-R", 4.1, tox=40.),
        con_full(3, 0, 2., "Drain"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=4,ivS=0.,ivE=10.,ivN=21,
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si logic-level power MOSFET 55V (IRLZ44N-like) planar VDMOS cell")

def T_rw_superjunction_600v():
    """Si superjunction MOSFET, 600 V class (CoolMOS-like), half cell.
    x=0 is the centre of an n pillar, x=Lx the centre of a p pillar. The
    drift is charge-balanced n and p pillars (Nd*Wn = Na*Wp, 3e15 x 3 um
    each): they deplete laterally, so the vertical field is nearly uniform and
    a 40 um drift holds 600 V with ~1/5 the on-resistance of a conventional
    design. The p-body sits on the p pillar; the n pillar reaches the surface
    (JFET region) where the gate-induced channel ends."""
    Lx,Ly,Lz = 6000., 44000., 400.; Nx,Ny,Nz = 44, 170, 4
    pill_d, body_d, src_d, sub_top = 41500., 1500., 300., 42000.
    body_x0, src_x0, src_x1 = 2000., 2400., 3400.
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,               "n",3e15, "n pillar"),
        Reg("Si",3000,0,0,Lx,pill_d,Lz,        "p",3e15, "p pillar (charge bal.)"),
        Reg("Si",body_x0,0,0,Lx,body_d,Lz,     "p",2e17, "p-body"),
        Reg("Si",0,sub_top,0,Lx,Ly,Lz,         "n",5e19, "n+ substrate"),
        Reg("Si",src_x0,0,0,src_x1,src_d,Lz,   "n",1e20, "n+ source"),
        Reg("Si",src_x1,0,0,Lx,src_d,Lz,       "p",1e19, "p+ body tap"),
    ]
    cons=[
        con_part(2, src_x0/Lx, 1.,          0.,1., 0, 0.,  "Source"),
        con_part(2, 0., (src_x0+200)/Lx,    0.,1., 2, 10., "Gate", 4.1, tox=80.),
        con_full(3, 0, 1., "Drain"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=2,ivS=0.,ivE=4.,ivN=17,
                dense_x=[body_x0,src_x0,src_x1,3000.],dense_y=[0.,src_d,body_d,sub_top],dense_z=[],
                desc="Si superjunction MOSFET 600V (CoolMOS-like) - charge-balanced pillars, half cell")

def T_rw_igbt_1200v():
    """Si field-stop trench IGBT, 1200 V / 40 A class (IKW40N120-like).
    n- base 1.5e14 cm^-3 x 105 um -> BV_ideal ~ 1245 V. The thin n buffer
    (field stop) lets the field terminate before reaching the collector,
    which is how modern 1200 V IGBTs stay much thinner than a non-punch-
    through design would need.  The collector is a lightly doped
    'transparent' p layer (3e17), as in field-stop IGBTs: it limits the hole
    injection efficiency, trading V_CE(sat) for switching loss.
    On-state current is several times a real part's (no inversion-layer
    mobility degradation, long lifetime); the plasma-dominated solve is the
    slowest in the library and stops converging near V_CE ~ 2 V."""
    Lx,Ly,Lz = 6000., 115000., 400.; Nx,Ny,Nz = 36, 200, 3
    body_d, src_d = 3000., 400.
    buf_top, coll_top = 111000., 114000.
    tw = 700.; tx0, tx1 = Lx/2.-tw, Lx/2.+tw; td = 5000.
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,            "n",1.5e14,"n- base (1200V)"),
        Reg("Si",0,buf_top,0,Lx,coll_top,Lz,"n",1e16,  "n buffer (field stop)"),
        Reg("Si",0,coll_top,0,Lx,Ly,Lz,     "p",3e17,  "p collector (transparent)"),
        Reg("Si",0,0,0,Lx,body_d,Lz,        "p",2e17,  "p-body"),
        Reg("Si",0,0,0,tx0,src_d,Lz,        "n",1e20,  "n+ emitter L"),
        Reg("Si",tx1,0,0,Lx,src_d,Lz,       "n",1e20,  "n+ emitter R"),
        Reg("Si",0,0,0,700,src_d,Lz,        "p",1e19,  "p+ tap L"),
        Reg("Si",Lx-700,0,0,Lx,src_d,Lz,    "p",1e19,  "p+ tap R"),
        Reg("SiO2",tx0,0,0,tx1,td,Lz,       "i",0.,    "trench oxide"),
    ]
    tox=100.
    cons=[
        con_part(2, 0., tx0/Lx,      0.,1., 0, 0.,  "Emitter-L"),
        con_part(2, tx1/Lx, 1.,      0.,1., 0, 0.,  "Emitter-R"),
        con_box(tx0+tox,0.,0.,tx1-tox,td-tox,Lz,(Lx,Ly,Lz), 2, 15., "Trench gate (n+ poly)", 4.1),
        con_full(3, 0, 1.5, "Collector"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=3,ivS=0.,ivE=1.8,ivN=19,
                dense_x=[tx0,tx0+100.,tx1-100.,tx1],dense_y=[src_d,body_d,td-100.,td],dense_z=[],
                desc="Si field-stop trench IGBT 1200V (IKW40N120-like) - oxide trench, buried gate (heavy: minutes)")

def T_rw_sic_mosfet_900v():
    """4H-SiC planar power MOSFET, 900 V / 65 mOhm class (C3M0065090D-like).
    Drift 1.8e16 cm^-3 x 8.5 um -> BV_ideal ~ 930 V. Note how much thinner
    and more heavily doped the drift is than a Si device of the same rating:
    that is the 10x critical-field advantage of SiC.
    Needs V_GS ~ 15-18 V because the SiC/oxide channel mobility is low."""
    Lx,Ly,Lz = 4000., 11500., 400.; Nx,Ny,Nz = 44, 130, 4
    body_d, src_d, sub_top = 800., 300., 9500.
    src_w, tap_w, jfet_hw  = 700., 300., 600.
    bx1_L = Lx/2.-jfet_hw; bx0_R = Lx/2.+jfet_hw
    regs=[
        Reg("4H-SiC",0,0,0,Lx,Ly,Lz,                     "n",1.8e16,"n drift (900V)"),
        Reg("4H-SiC",0,0,0,bx1_L,body_d,Lz,              "p",1e17,  "p-body L"),
        Reg("4H-SiC",bx0_R,0,0,Lx,body_d,Lz,             "p",1e17,  "p-body R"),
        Reg("4H-SiC",0,sub_top,0,Lx,Ly,Lz,               "n",1e19,  "n+ substrate"),
        Reg("4H-SiC",0,0,0,tap_w,src_d,Lz,               "p",1e19,  "p+ tap L"),
        Reg("4H-SiC",Lx-tap_w,0,0,Lx,src_d,Lz,           "p",1e19,  "p+ tap R"),
        Reg("4H-SiC",tap_w,0,0,tap_w+src_w,src_d,Lz,     "n",1e20,  "n+ source L"),
        Reg("4H-SiC",Lx-tap_w-src_w,0,0,Lx-tap_w,src_d,Lz,"n",1e20, "n+ source R"),
    ]
    sf=(src_w+tap_w)/Lx
    cons=[
        con_part(2, 0., sf,                              0.,1., 0, 0.,  "Source-L"),
        con_part(2, 1.-sf, 1.,                           0.,1., 0, 0.,  "Source-R"),
        con_part(2, (tap_w+src_w-200)/Lx, (bx1_L+200)/Lx,    0.,1., 2, 18., "Gate-L", 4.1, tox=50.),
        con_part(2, (bx0_R-200)/Lx, (Lx-tap_w-src_w+200)/Lx, 0.,1., 2, 18., "Gate-R", 4.1, tox=50.),
        con_full(3, 0, 2., "Drain"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=4,ivS=0.,ivE=10.,ivN=21,
                dense_x=[],dense_y=[],dense_z=[],
                desc="4H-SiC planar MOSFET 900V/65mOhm (C3M0065090D-like)")

def T_rw_sic_jbs_1200v():
    """4H-SiC Junction Barrier Schottky rectifier, 1200 V (C4D-class).
    Drift 1.2e16 cm^-3 x 12 um -> BV_ideal ~ 1400 V. The p+ grid pinches off
    the Schottky region under reverse bias, so leakage is far below a plain
    Schottky while forward drop stays Schottky-like (no minority injection,
    hence essentially zero reverse-recovery charge)."""
    Lx,Ly,Lz = 4000., 14500., 500.; Nx,Ny,Nz = 40, 130, 4
    sub_top, pp_d, pp_w = 12000., 500., 500.
    regs=[Reg("4H-SiC",0,0,0,Lx,Ly,Lz,          "n",1.2e16,"n drift (1200V)"),
          Reg("4H-SiC",0,sub_top,0,Lx,Ly,Lz,    "n",5e18,  "n+ substrate")]
    for cx in (700., 2000., 3300.):
        regs.append(Reg("4H-SiC",cx-pp_w/2,0,0,cx+pp_w/2,pp_d,Lz,
                        "p",1e18,f"p+ JBS grid @{cx:.0f}nm"))
    cons=[con_full(2,1,0.,"Schottky Anode",5.00), con_full(3,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-1400.,ivE=2.,ivN=33,bv_thr=0.05,
                dense_x=[],dense_y=[],dense_z=[],
                desc="4H-SiC JBS rectifier 1200V (C4D-class) - phi_B~1.2 eV")

def T_rw_gan_hemt_rf():
    """GaN RF power HEMT, 28 V / 10 W class (CGH40010-like), depletion mode.
    Lg=0.25 um, Lgd=1.55 um. 20 nm undoped Al0.25GaN barrier on GaN with the
    polarisation sheet charge (+1.2e13 q/cm^2) inducing the 2DEG.
    Ni Schottky gate, phi_m = 5.1 eV; pinch-off near -3...-4 V.
    No field plate here: at V_D = 28 V the peak field at the drain-side gate
    edge (~7 MV/cm) exceeds GaN's ~3.3 MV/cm critical field - real parts
    spread it with gate/source field plates (and there is no avalanche model
    to show the consequence)."""
    Lx,Ly,Lz = 3500., 600., 400.; Nx,Ny,Nz = 56, 54, 4
    bar=20.; src_w, drn_x = 500., 3000.
    regs=[
        Reg("GaN",0,0,0,Lx,Ly,Lz,              "p",1e15, "GaN buffer (SI)"),
        Reg("AlGaN",0,0,0,Lx,bar,Lz,           "n",1e10, "Al0.25GaN barrier"),
        Reg("GaN",0,0,0,src_w,80,Lz,           "n",5e19, "n+ source ohmic"),
        Reg("GaN",drn_x,0,0,Lx,80,Lz,          "n",5e19, "n+ drain ohmic"),
    ]
    cons=[
        con_part(2, 0., src_w/Lx,        0.,1., 0, 0.,  "Source"),
        con_part(2, drn_x/Lx, 1.,        0.,1., 0, 28., "Drain"),
        con_part(2, 1200/Lx, 1450/Lx,    0.,1., 1, -2., "Gate (Ni)", 5.1),
    ]
    sheets=[Sheet("y",bar,0.,Lx,0.,Lz,1.2e13,"AlGaN/GaN polarisation")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,sheets=sheets,
                ivc=1,ivS=0.,ivE=40.,ivN=21,
                dense_x=[src_w,1200.,1450.,drn_x],dense_y=[0.,bar,80.],dense_z=[],
                desc="GaN RF HEMT 28V/10W (CGH40010-like) - Lg=0.25um, polarisation 2DEG")
def T_rw_gan_hemt_650v():
    """GaN-on-Si e-mode power HEMT, 650 V class (GS66508-like).
    Enhancement mode by a p-GaN gate: 70 nm Mg-doped GaN (active Na ~1e18)
    plus the -sigma of the p-GaN/AlGaN interface lift the conduction band and
    deplete the 2DEG at V_G=0; elsewhere the surface is passivated (Al2O3
    stands in for SiN).  The gate contact to the p-GaN is modelled as OHMIC
    (GIT-like): with a Schottky metal the p-GaN floats between two back-to-back
    diodes and its DC potential is set by gate leakage, which is not modelled.
    Expect V_th ~ 1.5-2 V and a rising gate (pin-diode) current above ~3 V.
    12 nm AlGaN barrier with a net polarisation charge of 8e12 q/cm^2 (typical
    of the Al~0.2 barriers used under p-GaN gates; a thicker/higher-Al barrier
    makes the dipole across it so strong that the channel stays on at V_G=0).
    Blocking voltage is set laterally by Lgd = 15 um (no field plates here,
    and no impact-ionisation model - fields can be read, BV cannot)."""
    Lx,Ly,Lz = 20000., 1200., 500.; Nx,Ny,Nz = 70, 60, 4
    cap=70.; bar=12.; src_w, drn_x = 1500., 18000.; g0,g1=2500.,3000.
    regs=[
        Reg("GaN",0,0,0,Lx,Ly,Lz,               "p",1e15, "GaN buffer (SI)"),
        Reg("Al2O3",0,0,0,Lx,cap,Lz,            "i",0.,   "passivation (SiN stand-in)"),
        Reg("GaN",g0,0,0,g1,cap,Lz,             "p",1e18, "p-GaN gate (Mg, active)"),
        Reg("AlGaN",0,cap,0,Lx,cap+bar,Lz,      "n",1e10, "Al0.25GaN barrier"),
        Reg("GaN",0,0,0,src_w,cap+100,Lz,       "n",5e19, "n+ source ohmic"),
        Reg("GaN",drn_x,0,0,Lx,cap+100,Lz,      "n",5e19, "n+ drain ohmic"),
    ]
    cons=[
        con_part(2, 0., src_w/Lx,          0.,1., 0, 0.,   "Source"),
        con_part(2, drn_x/Lx, 1.,          0.,1., 0, 1., "Drain"),
        con_part(2, g0/Lx, g1/Lx,          0.,1., 0, 0.,   "Gate (p-GaN)"),
    ]
    # +sigma at AlGaN/GaN everywhere; under the gate also the -sigma of the
    # epitaxial p-GaN/AlGaN interface (in the access regions the AlGaN surface
    # -sigma is compensated by surface donors, the usual Ibbetson picture).
    # Together with the ionised Mg this lifts the channel above E_F: e-mode.
    sheets=[Sheet("y",cap+bar,0.,Lx,0.,Lz,8e12,"AlGaN/GaN polarisation"),
            Sheet("y",cap,g0,g1,0.,Lz,-8e12,"p-GaN/AlGaN polarisation")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,sheets=sheets,
                ivc=2,ivS=0.,ivE=3.,ivN=16,ivout="Drain",
                dense_x=[src_w,g0,g1,drn_x],dense_y=[cap,cap+bar,cap+100],dense_z=[],
                desc="GaN e-mode HEMT 650V (GS66508-like) - p-GaN gate (ohmic), polarisation 2DEG; sweep V_G")
def T_rw_ir_led_870():
    """AlGaAs/GaAs double-heterostructure IR emitter, ~870 nm (TSAL6200-like).
    GaAs active layer (Eg=1.42 eV -> 873 nm) confined between wider-gap
    Al(0.3)GaAs claddings, which is what gives a DH LED its high internal
    efficiency: carriers are trapped in the active layer instead of diffusing
    away. Used in every IR remote control and optocoupler."""
    Lx,Ly,Lz = 1000., 2300., 1000.; Nx,Ny,Nz = 8, 120, 6
    regs=[
        Reg("GaAs",  0,0,0,Lx,Ly,Lz,        "n",2e18,"n+ GaAs substrate"),
        Reg("AlGaAs",0,0,0,Lx,500,Lz,       "p",1e18,"p-AlGaAs cladding"),
        Reg("GaAs",  0,500,0,Lx,800,Lz,     "n",1e15,"GaAs active (873nm)"),
        Reg("AlGaAs",0,800,0,Lx,1300,Lz,    "n",1e18,"n-AlGaAs cladding"),
    ]
    cons=[con_full(2,0,0.,"Anode (p)"), con_full(3,0,0.,"Cathode (n)")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-3.,ivE=1.8,ivN=41,
                dense_x=[],dense_y=[],dense_z=[],
                desc="AlGaAs/GaAs DH IR LED ~870nm (TSAL6200-like)")

def T_rw_solar_perc():
    """Si PERC solar cell (thin-wafer variant).
    n+ emitter / p base / p+ rear BSF, with a front finger contact and a full
    rear contact. Base thinned to 60 um (real PERC wafers are 160-180 um) to
    keep the vertical aspect ratio tractable; the emitter and BSF fields,
    which set Voc, are unaffected.
    NOTE: there is no optical generation model, so this shows the dark
    electrostatics and dark I-V, not a photo I-V curve."""
    Lx,Ly,Lz = 3000., 61500., 500.; Nx,Ny,Nz = 20, 170, 4
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,           "p",1.5e16,"p base absorber"),
        Reg("Si",0,0,0,Lx,400,Lz,          "n",1e19,  "n+ front emitter"),
        Reg("Si",0,60500,0,Lx,Ly,Lz,       "p",5e18,  "p+ rear BSF"),
    ]
    cons=[con_part(2,0.35,0.65,0.,1.,0,0.,"Front finger (n+)"),
          con_full(3,0,0.,"Rear contact (p+)")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=-0.2,ivE=0.75,ivN=39,
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si PERC solar cell (thin wafer) - dark electrostatics only")

def T_rw_thyristor_800v():
    """Si thyristor / SCR, 800 V class (BT151-like).
    Four-layer p+ n- p n+ structure: two of the three junctions block in
    either polarity until a gate pulse turns on the regenerative npn/pnp pair.
    n- base 3e14 cm^-3 x 60 um -> BV_ideal ~ 740 V.
    Default: forward biased at V_AK = 2 V, gate at 0 V (blocking: J2 holds the
    voltage, the n- base floats).  The sweep raises the gate voltage: J3
    injects, the npn and pnp gains add up and the device LATCHES near
    V_G ~ 0.65 V (the anode current jumps by ~3 decades; with a voltage-
    forced gate the latched p-base then sinks current through the gate).
    For the 400 V blocking state set Anode = 400 V (slower: the floating
    n- base is resolved by a separate balance step)."""
    Lx,Ly,Lz = 10000., 72000., 500.; Nx,Ny,Nz = 40, 170, 4
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,           "n",3e14,"n- base (800V)"),
        Reg("Si",0,0,0,Lx,10000,Lz,        "p",1e16,"p base (to surface around cathode)"),
        Reg("Si",0,0,0,6000,2000,Lz,       "n",1e19,"n+ cathode"),
        Reg("Si",8000,0,0,Lx,2000,Lz,      "p",1e19,"p+ gate contact"),
        Reg("Si",0,70000,0,Lx,Ly,Lz,       "p",1e19,"p+ anode"),
    ]
    cons=[
        con_part(2, 0., 6000/Lx,     0.,1., 0, 0.,   "Cathode"),
        con_part(2, 8000/Lx, 1.,     0.,1., 0, 0.,   "Gate"),
        con_full(3, 0, 2., "Anode"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=0.,ivE=1.0,ivN=21,ivout="Anode",
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si thyristor/SCR 800V (BT151-like) - 4-layer pnpn; gate sweep shows latching")

# ══════════════════════════════════════════════════════════════
#  CMOS LOGIC
# ══════════════════════════════════════════════════════════════
def T_cmos_inverter():
    """CMOS inverter, 0.25 um class (VDD = 2.5 V), 2-D cross-section.

    Twin wells (p-well | n-well, 5e17, 800 nm deep) in a 5e16 p-substrate,
    separated by a 350 nm shallow-trench isolation (STI).  NMOS: n+ poly gate
    (phi_m 4.05 eV); PMOS: p+ poly gate (5.17 eV); both on 5 nm SiO2, L = 250 nm,
    abrupt 150 nm deep n+/p+ source/drain.  Long-channel estimate
    V_th = V_FB + 2 phi_F + Q_dep/C_ox = +0.47 V (NMOS) / -0.47 V (PMOS);
    the simulated values are +0.54 V / -0.54 V (+0.48 / -0.49 V with the
    classical models only: quantum confinement and surface mobility shift the
    max-g_m threshold).

    Wiring (nets): both gates = "Vin", both drains = "Out", the NMOS source +
    p+ well tap and the substrate contact = "GND", the PMOS source + n+ well
    tap = "VDD".

    The sweep is the VOLTAGE TRANSFER CURVE: Vin 0 -> 2.5 V with the Out net
    FLOATING - at every point Out is solved for zero current (an unloaded
    output) and the steep part is refined automatically.  The I-V tab shows
    Vout(Vin), the gain (~21), the switching threshold V_M (~1.06 V), the
    noise margins and the short-circuit supply-current peak (~15 uA per cell).

    Both transistors have the same width (one cross-section, 200 nm deep in
    z), so the pull-up is weaker (surface mobility mu_p ~ mu_n/2.5) and V_M
    sits below VDD/2; real cells make W_p ~ 2-3 W_n to centre it.  Channel
    mobility includes surface (Lombardi) degradation and the inversion
    layers are quantum-corrected (MLDA).  A sweep takes 1-3 minutes (each
    point solves the output node iteratively)."""
    Lx,Ly,Lz = 3000., 1300., 200.; Nx,Ny,Nz = 96, 44, 3
    sd, sti_d, well, xw = 150., 350., 800., 1500.
    tapN,srcN,gN,drnN = (0.,200.),(200.,450.),(450.,700.),(700.,1300.)
    sti = (1300.,1700.)
    drnP,gP,srcP,tapP = (1700.,2300.),(2300.,2550.),(2550.,2800.),(2800.,3000.)
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,              "p",5e16,"p substrate"),
        Reg("Si",0,0,0,xw,well,Lz,            "p",5e17,"p-well"),
        Reg("Si",xw,0,0,Lx,well,Lz,           "n",5e17,"n-well"),
        Reg("Si",tapN[0],0,0,tapN[1],sd,Lz,   "p",1e20,"p+ well tap"),
        Reg("Si",srcN[0],0,0,srcN[1],sd,Lz,   "n",1e20,"n+ source (NMOS)"),
        Reg("Si",drnN[0],0,0,drnN[1],sd,Lz,   "n",1e20,"n+ drain (NMOS)"),
        Reg("Si",drnP[0],0,0,drnP[1],sd,Lz,   "p",1e20,"p+ drain (PMOS)"),
        Reg("Si",srcP[0],0,0,srcP[1],sd,Lz,   "p",1e20,"p+ source (PMOS)"),
        Reg("Si",tapP[0],0,0,tapP[1],sd,Lz,   "n",1e20,"n+ well tap"),
        Reg("SiO2",sti[0],0,0,sti[1],sti_d,Lz,"i",0.,  "shallow trench isolation"),
    ]
    fx=lambda a,b:(a/Lx,b/Lx)
    VDD=2.5
    cons=[
        con_part(2,*fx(tapN[0],srcN[1]),0.,1.,0,0.,  "GND",                    net="GND"),
        con_part(2,*fx(*gN),0.,1.,2,0.,  "Gate N (n+ poly)",4.05,tox=5.,       net="Vin"),
        con_part(2,*fx(*drnN),0.,1.,0,VDD,"Out N",                             net="Out"),
        con_part(2,*fx(*drnP),0.,1.,0,VDD,"Out P",                             net="Out"),
        con_part(2,*fx(*gP),0.,1.,2,0.,  "Gate P (p+ poly)",5.17,tox=5.,       net="Vin"),
        con_part(2,*fx(srcP[0],tapP[1]),0.,1.,0,VDD,"VDD"),
        con_full(3,0,0.,"Substrate",                                            net="GND"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=0.,ivE=VDD,ivN=26,ivnet="Vin",ivfloat="Out",ivout="VDD",models=63,
                dense_x=[tapN[1],srcN[1],gN[1],drnN[1],xw,drnP[1],gP[1],srcP[1]],
                dense_y=[0.,sd,sti_d,well],dense_z=[],
                desc="CMOS inverter 0.25um (VDD 2.5V) - twin well + STI; transfer curve with floating output")

def T_finfet():
    """Bulk n-FinFET, 22 nm node class (tri-gate), full 3-D, high-k/metal gate.

    x = source -> drain along the fin, y = depth (fin top at y = 5 nm),
    z = across the fin.  Fin 10 nm wide and 35 nm tall above the STI,
    L_g = 30 nm.  Gate stack on the top and both sidewalls: 0.5 nm SiO2
    interfacial layer + 2.5 nm HfO2 (k = 22), EOT = 0.94 nm, and a TiN-like
    metal gate (phi_m 4.50 eV) - three box electrodes on the net "Gate",
    swept as one terminal.  The stack is meshed node by node (the template's
    "nodes"/"clear" keys), so its layers have exactly these thicknesses.

    Physics on (Solver page): quantum confinement (MLDA - the inversion layer
    sits ~1 nm off the oxide, which raises V_th and thins the effective
    fin), Lombardi surface mobility and gate tunnelling through the stack.

    The fin is nearly undoped (1e16) - the gate, not the channel doping, sets
    V_th - and a 3e18 punch-through stopper in the sub-fin (below the gate's
    reach, 2 nm under the S/D) blocks the leakage path under the fin.  n+
    (1e20) source/drain fins with contacts on their end faces; 1e18
    p-substrate, body contact at the bottom.

    Default V_D = 0.8 V; the sweep is the transfer curve I_D(V_G), 0 -> 0.8 V.
    Simulated: SS = 67 mV/dec, V_th = 0.29 V (constant current 100 nA x W/L,
    W_eff = 80 nm), I_on ~ 48 uA per fin, I_on/I_off ~ 3e6, DIBL ~ 31 mV/V
    (set Drain = 0.05 V and sweep again).  Gate leakage through the HfO2
    stack ~3e-12 A per fin at V_G = V_D = 0.8 V (~0.1 A/cm^2); a SiO2 gate of
    the same EOT would leak ~1000x more (see the gate-leakage templates).
    The Structure tab opens on the YZ section through the gate.

    Not modelled: ballistic/quasi-ballistic transport and strain, so I_on is
    a drift-diffusion (velocity-saturated) estimate."""
    Lx,Ly,Lz = 90., 110., 40.; Nx,Ny,Nz = 30, 30, 16
    g0,g1 = 30., 60.                 # gate: x in [30,60) -> L_g = 30 nm
    zf0,zf1 = 15., 25.05             # fin: Si nodes z = 15.0 ... 25.0
    ytop, ysti, ypts, ybase = 5., 40., 42., 70.
    pm = 4.50
    # stack interfaces (node-pair midpoints): Si|IL at 4.95 / 14.95 / 25.05,
    # IL|HfO2 at 4.45 / 14.45 / 25.55, metal surface (metal node) at 1.95 /
    # 11.95 / 28.05  ->  IL 0.5 nm, HfO2 2.5 nm on all three sides
    yi, ym = 4.45, 1.95
    zl_i, zl_m = 14.45, 11.95
    zr_i, zr_m = 25.55, 28.05
    regs=[
        Reg("SiO2",0,0,0,Lx,Ly,Lz,           "i",0.,  "oxide: spacer, STI"),
        Reg("Si",0,ybase,0,Lx,Ly,Lz,         "p",1e18,"p substrate"),
        Reg("Si",0,ytop,zf0,Lx,ybase,zf1,    "p",1e16,"fin (undoped)"),
        Reg("Si",0,ypts,zf0,Lx,ybase,zf1,    "p",3e18,"punch-through stopper"),
        Reg("Si",0,ytop,zf0,g0,ysti,zf1,     "n",1e20,"n+ source fin"),
        Reg("Si",g1,ytop,zf0,Lx,ysti,zf1,    "n",1e20,"n+ drain fin"),
        # gate stack under the gate: HfO2 wrapping the fin, SiO2 interfacial layer
        # (the metal-surface nodes belong to the HfO2 regions: that edge is HfO2)
        Reg("HfO2",g0,ym,zl_m,g1,yi,zr_m+0.01,            "i",0.,"HfO2 (top)"),
        Reg("HfO2",g0,ym,zl_m,g1,ysti,zl_i,               "i",0.,"HfO2 (left)"),
        Reg("HfO2",g0,ym,zr_i,g1,ysti,zr_m+0.01,          "i",0.,"HfO2 (right)"),
        Reg("SiO2",g0,yi,zl_i,g1,ytop,zr_i,               "i",0.,"interfacial SiO2 (top)"),
        Reg("SiO2",g0,yi,zl_i,g1,ysti,zf0,                "i",0.,"interfacial SiO2 (left)"),
        Reg("SiO2",g0,yi,zf1,g1,ysti,zr_i,                "i",0.,"interfacial SiO2 (right)"),
    ]
    L=(Lx,Ly,Lz)
    # X-faces: i -> Y, j -> Z (fractions): contacts on the S/D fin end faces
    cons=[
        Con(0, ytop/Ly, ysti/Ly, zf0/Lz, zf1/Lz, 0, 0.,  4.05, "Source"),
        Con(1, ytop/Ly, ysti/Ly, zf0/Lz, zf1/Lz, 0, 0.8, 4.05, "Drain"),
        con_box(g0,0.,0.,   g1,ym+0.01,Lz,   L, 2, 0.8, "Gate top",   pm, net="Gate"),
        con_box(g0,0.,0.,   g1,ysti,zl_m+0.01,L, 2, 0.8, "Gate left",  pm, net="Gate"),
        con_box(g0,0.,zr_m-0.01,g1,ysti,Lz,  L, 2, 0.8, "Gate right", pm, net="Gate"),
        con_full(3,0,0.,"Body"),
    ]
    nodes=dict(y=[ym,3.2,4.4,4.5,4.9,5.0],
               z=[zl_m,13.2,14.4,14.5,14.9,15.0,25.0,25.1,25.5,25.6,26.8,zr_m])
    clear=dict(y=[(ym,5.0)], z=[(zl_m,15.0),(25.0,zr_m)])
    refine=[("y",ytop,+1),("z",15.,+1),("z",25.,-1)]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=2,ivS=0.,ivE=0.8,ivN=17,ivnet="Gate",ivout="Drain",models=63,
                nodes=nodes,clear=clear,refine=refine,view=("YZ",45.),wl=80./30.,cut="X",extrude="x",
                dense_x=[g0,g1],dense_y=[ytop,ysti,ypts,ybase],dense_z=[15.,25.],
                desc="Bulk n-FinFET 22nm class - tri-gate, 10nm fin, Lg 30nm, SiO2/HfO2 stack EOT 0.94nm; I_D(V_G) at V_D=0.8V")

# ══════════════════════════════════════════════════════════════
#  GATE DIELECTRICS: tunnelling leakage
# ══════════════════════════════════════════════════════════════
def _nmos_leak(stack, pm, Na, ivS, ivE, ivN, desc, tox_note=""):
    """Planar n-MOSFET cross-section for gate-current studies: n+ source and
    drain (grounded) supply the inversion layer, the p-body contact is at the
    bottom and the gate is an oxide (Robin) gate on the top face with the
    given layer stack."""
    Lx,Ly,Lz = 300., 150., 20.; Nx,Ny,Nz = 36, 30, 3
    regs=[Reg("Si",0,0,0,Lx,Ly,Lz,       "p",Na,  "p-body"),
          Reg("Si",0,0,0,100,50,Lz,      "n",1e20,"n+ source"),
          Reg("Si",200,0,0,Lx,50,Lz,     "n",1e20,"n+ drain")]
    cons=[con_part(2,0.,0.30,0.,1.,0,0.,"Source"),
          con_part(2,0.70,1.,0.,1.,0,0.,"Drain"),
          con_part(2,100./Lx,200./Lx,0.,1.,2,1.0,"Gate",pm,stack=stack),   # no S/D overlap
          con_full(3,0,0.,"Body")]
    c=cons[2]; c.tox=round(c.eot(),4)
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=2,ivS=ivS,ivE=ivE,ivN=ivN,ivnet="Gate",ivout="Gate",models=63,
                dense_x=[100.,200.],dense_y=[0.,50.],dense_z=[],desc=desc)

def T_finfet_rounded():
    """Bulk n-FinFET with a REALISTIC fin profile: tapered sidewalls and
    rounded top corners (22/14 nm node class), full 3-D, high-k/metal gate.

    Same device as "FinFET (bulk tri-gate)" except the fin shape: 8 nm wide
    at the top, 12 nm at the STI surface (sidewall angle 86.7 deg, as etched
    fins are), top corners rounded with r = 3 nm.  The 0.5 nm SiO2 + 2.5 nm
    HfO2 stack is conformal (constant thickness all round, corners rounded
    r + t) and the metal gate (phi_m 4.50 eV) fills everything outside it.
    L_g = 30 nm, n+ source/drain, 3e18 punch-through stopper, body contact
    at the bottom; V_D = 0.8 V, sweep I_D(V_G) 0 -> 0.8 V.

    The shapes need the TRIANGULAR (prism) mesh (the template opens with
    it): sloped and rounded interfaces are meshed as they are, with node
    pairs straddling each interface and the stack layers aligned along the
    normals.  Rounded corners weaken the corner effect of a sharp fin (field
    crowding that turns the corners on first) and the taper widens the fin
    towards its base, so compare V_th, SS and I_on with the sharp template.
    On the rectangular mesh the same outlines become staircases."""
    Lx,Ly,Lz = 90., 110., 40.; Nx,Ny,Nz = 30, 40, 26
    g0,g1 = 30., 60.
    ytop, ysti, ypts, ybase = 5., 40., 42., 70.
    zc, wtop, wsti, rc = 20., 8., 12., 3.
    slope = 0.5*(wsti-wtop)/(ysti-ytop)
    hb = 0.5*wsti+slope*(ybase-ysti)
    fin = SM.Shape("x", [(ytop,zc-0.5*wtop),(ytop,zc+0.5*wtop),(ybase,zc+hb),(ybase,zc-hb)], [rc,rc,0.,0.])
    il, hk = fin.offset(0.5), fin.offset(3.0)
    regs=[
        Reg("SiO2",0,0,0,Lx,Ly,Lz,           "i",0.,  "oxide: spacer, STI"),
        Reg("Si",0,ybase,0,Lx,Ly,Lz,         "p",1e18,"p substrate"),
        Reg("HfO2",g0,0,0,g1,ysti,Lz,        "i",0.,  "HfO2",            shape=hk),
        Reg("SiO2",g0,0,0,g1,ysti,Lz,        "i",0.,  "interfacial SiO2",shape=il),
        Reg("Si",0,ytop,0,Lx,ybase,Lz,       "p",1e16,"fin (undoped)",   shape=fin),
        Reg("Si",0,ypts,0,Lx,ybase,Lz,       "p",3e18,"punch-through stopper",shape=fin),
        Reg("Si",0,ytop,0,g0,ysti,Lz,        "n",1e20,"n+ source fin",   shape=fin),
        Reg("Si",g1,ytop,0,Lx,ysti,Lz,       "n",1e20,"n+ drain fin",    shape=fin),
    ]
    L=(Lx,Ly,Lz)
    zs0,zs1 = zc-0.5*wsti, zc+0.5*wsti
    cons=[
        Con(0, ytop/Ly, ysti/Ly, zs0/Lz, zs1/Lz, 0, 0.,  4.05, "Source"),
        Con(1, ytop/Ly, ysti/Ly, zs0/Lz, zs1/Lz, 0, 0.8, 4.05, "Drain"),
        con_box(g0,0.,0., g1,ysti,Lz, L, 2, 0.8, "Gate", 4.50, shape=hk.outside()),
        con_full(3,0,0.,"Body"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=2,ivS=0.,ivE=0.8,ivN=17,ivnet="Gate",ivout="Drain",models=63,
                mesh="tri",extrude="x",view=("YZ",45.),wl=80./30.,cut="X",
                dense_x=[g0,g1],dense_y=[ytop,ysti,ypts,ybase],dense_z=[zs0,zc,zs1],
                desc="Bulk n-FinFET with tapered, rounded fin (8/12nm, r 3nm), conformal SiO2/HfO2 stack; I_D(V_G) at V_D=0.8V")

def T_gaa_nanowire():
    """Gate-all-around Si nanowire n-FET (nanosheet-era class), full 3-D.

    x = source -> drain along the wire; y, z = its cross-section.  A round Si
    wire 10 nm in diameter, L_g = 20 nm, wrapped completely by 0.6 nm SiO2
    + 1.9 nm HfO2 (EOT 0.94 nm) and a TiN-like metal gate (phi_m 4.55 eV).
    n+ (1e20) source/drain wire sections abrupt at the gate edges, 1e15 p
    channel, SiO2 spacers around the wire outside the gate, contacts on the
    wire end faces.  V_D = 0.7 V; the sweep is I_D(V_G) 0 -> 0.7 V.

    This is the geometry the TRIANGULAR (prism) mesh is for (the template
    opens with it): the wire, both stack layers and the gate surface are
    meshed as true circles - a node pair straddles every interface and the
    stack's nodes line up along each radius, so the tunnelling paths run
    straight through it.  Switch the Device page to the rectangular mesh to
    see the staircase version of the same wire.

    With a gate all around, the channel electrostatics are better than in a
    tri-gate fin of similar size: SS close to 60 mV/dec and a small DIBL
    (set Drain = 0.05 V and sweep again).  Quantum confinement (MLDA) moves
    the electrons ~1 nm off the whole oxide ring; Lombardi surface mobility
    and gate tunnelling are on.  Not modelled: ballistic transport, strain
    and the wire's 1-D subband quantisation beyond the MLDA density
    correction."""
    Lx,Ly,Lz = 80., 22., 22.; Nx,Ny,Nz = 30, 34, 34
    g0,g1 = 30., 50.
    uc = vc = 11.; R = 5.
    wire = SM.Shape.circle("x", uc, vc, R)
    il, hk = wire.offset(0.6), wire.offset(2.5)
    regs=[
        Reg("SiO2",0,0,0,Lx,Ly,Lz,       "i",0.,  "spacer oxide"),
        Reg("HfO2",g0,0,0,g1,Ly,Lz,      "i",0.,  "HfO2",            shape=hk),
        Reg("SiO2",g0,0,0,g1,Ly,Lz,      "i",0.,  "interfacial SiO2",shape=il),
        Reg("Si",0,0,0,Lx,Ly,Lz,         "p",1e15,"channel (wire)",  shape=wire),
        Reg("Si",0,0,0,g0,Ly,Lz,         "n",1e20,"n+ source",       shape=wire),
        Reg("Si",g1,0,0,Lx,Ly,Lz,        "n",1e20,"n+ drain",        shape=wire),
    ]
    L=(Lx,Ly,Lz)
    f0,f1 = (uc-R)/Ly, (uc+R)/Ly
    cons=[
        Con(0, f0, f1, f0, f1, 0, 0.,  4.05, "Source"),
        Con(1, f0, f1, f0, f1, 0, 0.7, 4.05, "Drain"),
        con_box(g0,0.,0., g1,Ly,Lz, L, 2, 0.7, "Gate", 4.55, shape=hk.outside()),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=2,ivS=0.,ivE=0.7,ivN=15,ivnet="Gate",ivout="Drain",models=63,
                mesh="tri",extrude="x",view=("YZ",40.),wl=np.pi*2*R/Lg if (Lg:=g1-g0) else 1.,cut="X",
                dense_x=[g0,g1],dense_y=[],dense_z=[],
                desc="GAA Si nanowire n-FET - 10nm round wire, Lg 20nm, SiO2/HfO2 stack wrapped all around; I_D(V_G) at V_D=0.7V")

def T_leak_sio2():
    """Gate leakage of a 1.2 nm SiO2 gate (planar n-MOSFET, 90 nm node class).

    Metal gate (phi_m 4.15 eV) on 1.2 nm SiO2 over a 1e18 p-body, n+ source
    and drain grounded; the sweep is V_G from -1.5 to +1.5 V and the plotted
    output is the GATE current (I-V tab: linear and semi-log, J_G at +-0.5 and
    +-1 V in the parameter list).

    At this thickness electrons tunnel directly through the oxide: at
    V_G = +1 V (inversion) J_G ~ 80 A/cm^2 (textbook range 1e2-1e3 A/cm^2 for
    1.2 nm SiO2), rising ~10x for every 0.2 nm less oxide.  The electron
    current comes from the inversion layer (source and drain feed it - see
    their currents); at negative V_G (accumulation, V_FB = -0.94 V) electrons
    tunnel from the gate into the substrate conduction band and holes from
    the accumulation layer into the gate (~30 A/cm^2 at -1 V).  The gate does
    not overlap source and drain (an overlap adds edge tunnelling, which
    dominates at negative V_G).  Compare the high-k template of the same EOT.

    Tunnelling model: Tsu-Esaki with WKB transmission through the stack
    (SiO2 barrier 3.1 eV, tunnelling mass 0.42 m0), supply from the
    inversion/accumulation layer's pressure on the interface; direct and
    Fowler-Nordheim tunnelling in one formula.  MLDA and Lombardi are on."""
    return _nmos_leak([("SiO2",1.2)],4.15,1e18,-1.5,1.5,31,
                      "NMOS gate leakage - 1.2nm SiO2 (direct tunnelling), sweep V_G, plot I_G")

def T_leak_hfo2():
    """Gate leakage of a high-k stack of the SAME EOT (1.2 nm): 0.5 nm SiO2
    interfacial layer + 3.95 nm HfO2 (k = 22), metal gate (phi_m 4.15 eV).

    Identical transistor to "Gate leakage: 1.2nm SiO2 NMOS" - same gate
    capacitance, same V_th and inversion charge - but the physical dielectric
    is 4.45 nm thick.  HfO2's conduction-band offset to Si is only 1.5 eV
    (SiO2: 3.1 eV) and its tunnelling mass is lower (0.18 m0), yet the
    extra thickness wins: J_G(+1 V) ~ 1.5e-2 A/cm^2, ~5000x below the SiO2
    gate; the 0.5 nm interfacial SiO2 now carries most of the tunnelling
    resistance.  This is direct tunnelling only: real stacks also leak through
    traps in the HfO2 (not modelled) and measured reductions are 1e2-1e4x.
    This is why the 45 nm node and all later ones moved to high-k/metal gates.

    Change the stack on the Contacts page (e.g. 'SiO2 0.8, HfO2 2.2') to see
    how the interfacial layer dominates the leakage."""
    return _nmos_leak([("SiO2",0.5),("HfO2",3.95)],4.15,1e18,-1.5,1.5,31,
                      "NMOS gate leakage - SiO2/HfO2 high-k stack, EOT 1.2nm; sweep V_G, plot I_G")

def T_fowler_nordheim():
    """Fowler-Nordheim tunnelling through 8 nm SiO2 (thick-oxide MOSFET).

    Direct tunnelling through 8 nm is negligible; once the oxide drop exceeds
    the 3.1 eV barrier the barrier turns triangular and electrons from the
    inversion layer tunnel into the oxide conduction band (Fowler-Nordheim):
    J = A E^2 exp(-B/E) with B ~ 240 MV/cm for SiO2 (m_ox = 0.42 m0).  The
    sweep takes V_G from 4 to 10 V (oxide field up to ~12 MV/cm, around the
    breakdown field of real SiO2): the gate current rises over ~8 decades,
    ~1e-9 A/cm^2 at 5 V to ~0.4 A/cm^2 at 10 V - the regime used to program
    flash memory cells, and the reason oxide fields above ~5 MV/cm degrade
    MOS gates in operation.  Plot ln(J/E^2) against 1/E (export the sweep)
    for the straight Fowler-Nordheim line."""
    return _nmos_leak([("SiO2",8.0)],4.15,5e17,4.,10.,25,
                      "Fowler-Nordheim tunnelling through 8nm SiO2 - sweep V_G 4-10V, plot I_G")

TEMPLATES = {
    # ── textbook / idealised structures ──────────────────────────
    "PN Diode":T_pn,"PIN Diode":T_pin,"Schottky":T_schottky,
    "NPN BJT":T_npn,"PNP BJT":T_pnp,"NMOS":T_nmos,"PMOS":T_pmos,
    "VDMOS":T_vdmos,"LDMOS":T_ldmos,"IGBT":T_igbt,
    "GaAs LED":T_led_gaas,
    "GaAs HBT":T_hbt,"Si Solar":T_solar,"4H-SiC PiN":T_sic_pn,
    "GaN HEMT":T_gan_hemt,"n-JFET":T_jfet,
    "SiC Trench MOS":T_trench_sic_mos,"Si Trench MOS":T_trench_si_mos,
    "SiC JBS Diode":T_sic_jbs,
    # ── real-world parts: diodes & rectifiers ───────────────────
    "RW: Si diode 100V (1N4148)":      T_rw_diode_1n4148,
    "RW: Si Schottky 40V (1N5819)":    T_rw_schottky_1n5819,
    "RW: Si Zener 5.1V (1N4733A)":     T_rw_zener_5v1,
    "RW: SiC JBS 1200V (C4D)":         T_rw_sic_jbs_1200v,
    # ── real-world parts: optoelectronics ───────────────────────
    "RW: Si PIN photodiode (BPW34)":   T_rw_pin_photodiode,
    "RW: Ge photodiode 1550nm":        T_rw_ge_photodiode,
    "RW: AlGaAs IR LED 870nm":         T_rw_ir_led_870,
    "RW: Si PERC solar cell":          T_rw_solar_perc,
    # ── real-world parts: transistors & switches ────────────────
    "RW: Si NPN BJT 45V (BC547)":      T_rw_bjt_bc547,
    "RW: Si MOSFET 55V (IRLZ44N)":     T_rw_mosfet_irlz44n,
    "RW: Si superjunction 600V":       T_rw_superjunction_600v,
    "RW: Si IGBT 1200V (IKW40N120)":   T_rw_igbt_1200v,
    "RW: SiC MOSFET 900V (C3M)":       T_rw_sic_mosfet_900v,
    "RW: GaN RF HEMT 28V (CGH40010)":  T_rw_gan_hemt_rf,
    "RW: GaN HEMT 650V (GS66508)":     T_rw_gan_hemt_650v,
    "RW: Si thyristor 800V (BT151)":   T_rw_thyristor_800v,
    # ── CMOS logic ───────────────────────────────────────────────
    "CMOS Inverter (0.25um)":          T_cmos_inverter,
    "FinFET (bulk tri-gate)":          T_finfet,
    "FinFET (tapered, rounded fin)":   T_finfet_rounded,
    "GAA nanowire FET (round)":        T_gaa_nanowire,
    # ── gate dielectrics ─────────────────────────────────────────
    "Gate leakage: 1.2nm SiO2 NMOS":   T_leak_sio2,
    "Gate leakage: HfO2 high-k NMOS":  T_leak_hfo2,
    "Fowler-Nordheim: 8nm SiO2 MOS":   T_fowler_nordheim,
}
# Menu categories (every template appears once; anything not listed here
# lands in "Other")
TEMPLATE_GROUPS = [
    ("Diodes & rectifiers", ["PN Diode","PIN Diode","Schottky","4H-SiC PiN","SiC JBS Diode",
        "RW: Si diode 100V (1N4148)","RW: Si Schottky 40V (1N5819)",
        "RW: Si Zener 5.1V (1N4733A)","RW: SiC JBS 1200V (C4D)"]),
    ("MOSFETs, CMOS & FinFET", ["NMOS","PMOS","CMOS Inverter (0.25um)","FinFET (bulk tri-gate)",
        "FinFET (tapered, rounded fin)","GAA nanowire FET (round)",
        "n-JFET"]),
    ("Gate dielectrics & tunnelling", ["Gate leakage: 1.2nm SiO2 NMOS","Gate leakage: HfO2 high-k NMOS",
        "Fowler-Nordheim: 8nm SiO2 MOS"]),
    ("Power switches", ["VDMOS","LDMOS","Si Trench MOS","SiC Trench MOS","IGBT",
        "RW: Si MOSFET 55V (IRLZ44N)","RW: Si superjunction 600V",
        "RW: Si IGBT 1200V (IKW40N120)","RW: SiC MOSFET 900V (C3M)"]),
    ("Bipolar & thyristors", ["NPN BJT","PNP BJT","GaAs HBT","RW: Si NPN BJT 45V (BC547)",
        "RW: Si thyristor 800V (BT151)"]),
    ("GaN HEMTs", ["GaN HEMT","RW: GaN RF HEMT 28V (CGH40010)","RW: GaN HEMT 650V (GS66508)"]),
    ("Optoelectronics & solar", ["GaAs LED","Si Solar","RW: Si PIN photodiode (BPW34)",
        "RW: Ge photodiode 1550nm","RW: AlGaAs IR LED 870nm","RW: Si PERC solar cell"]),
]
def template_groups():
    """[(category, [names])] covering every template exactly once."""
    seen=set(); out=[]
    for cat,names in TEMPLATE_GROUPS:
        ns=[n for n in names if n in TEMPLATES and n not in seen]; seen.update(ns)
        if ns: out.append((cat,ns))
    rest=[n for n in TEMPLATES if n not in seen]
    if rest: out.append(("Other",rest))
    return out
def template_doc(name):
    """Long description of a template: its function's docstring (or desc)."""
    fn=TEMPLATES.get(name)
    doc=(fn.__doc__ or "").strip() if fn else ""
    if not doc:
        try: doc=fn()["desc"]
        except Exception: doc=""
    import textwrap
    lines=doc.splitlines()
    if len(lines)>1:
        doc=lines[0].strip()+"\n"+textwrap.dedent("\n".join(lines[1:]))
    # reflow: join wrapped lines of a paragraph; keep blank lines, indented
    # lines and list items on their own line
    out=[]; prev_ind=False
    for ln in doc.split("\n"):
        st=ln.strip(); ind=ln.startswith((" ","\t"))
        if not st: out.append(""); prev_ind=False; continue
        item=st[:2] in ("- ","* ")
        if out and out[-1] and not item and ((not ind and not prev_ind) or
                (ind and prev_ind and not out[-1].rstrip().endswith((".",":")))):
            out[-1]+=" "+st
        else: out.append(("   " if ind else "")+st)
        prev_ind=ind
    return "\n".join(out)

# ══════════════════════════════════════════════════════════════
#  DEVICE BUILDER (shared by the GUI and the validation scripts)
# ══════════════════════════════════════════════════════════════
def insert_layer(xs, pos_m, h0=0.5e-9, growth=1.3, side=0):
    """Insert a geometric boundary layer of nodes at pos (m): spacings h0,
    h0*g, h0*g^2 ... away from pos (side=+1 towards larger x, -1 towards
    smaller, 0 both) until the existing grid is at least as fine. Resolves
    inversion layers / 2DEGs (a few nm) that the 8x graded grid cannot."""
    arr=np.asarray(xs,dtype=float); L0,L1=arr[0],arr[-1]
    pts=set(arr.tolist())
    if L0<pos_m<L1: pts.add(pos_m)
    for sg in ([1,-1] if side==0 else [side]):
        x=pos_m; dd=h0
        while True:
            x=x+sg*dd
            if not (L0<x<L1): break
            near=np.min(np.abs(arr-x))
            if near<0.7*dd: break           # existing grid already this fine
            pts.add(x); dd*=growth
    out=np.array(sorted(pts))
    keep=np.concatenate(([True],np.diff(out)>1e-13))
    return out[keep]

def make_grid(t, graded=True):
    """Node coordinates (m) for a template/device dict."""
    Lx,Ly,Lz=t["Lx"],t["Ly"],t["Lz"]; Nx,Ny,Nz=t["Nx"],t["Ny"],t["Nz"]
    if graded:
        # Auto-detected junction/material planes merged with manual dense lists.
        # An axis without junctions stays uniform (no spurious centre refinement).
        ax,ay,az=detect_junctions(t["regs"],Lx,Ly,Lz)
        dx=sorted(set(list(t.get("dense_x",[]))+ax))
        dy=sorted(set(list(t.get("dense_y",[]))+ay))
        dz=sorted(set(list(t.get("dense_z",[]))+az))
        xs=graded_grid(Lx,Nx,dx) if dx else uniform_grid(Lx,Nx)
        ys=graded_grid(Ly,Ny,dy) if dy else uniform_grid(Ly,Ny)
        zs=graded_grid(Lz,Nz,dz) if dz else uniform_grid(Lz,Nz)
    else:
        xs=uniform_grid(Lx,Nx); ys=uniform_grid(Ly,Ny); zs=uniform_grid(Lz,Nz)
    if graded:
        g=[xs,ys,zs]; Ls=[Lx,Ly,Lz]
        # oxide-gated faces: inversion/accumulation layer at the surface
        for c in t.get("cons",[]):
            if c.bc==2 and c.face<=5:
                ax=c.face//2; pos=0. if c.face%2==0 else Ls[ax]*1e-9
                g[ax]=insert_layer(g[ax],pos,0.5e-9,1.3,+1 if c.face%2==0 else -1)
        # oxide/semiconductor interfaces (trench sidewalls, buried oxides):
        # refine on the semiconductor side of every insulator face
        for r in t.get("regs",[]):
            if r.mat not in INSULATORS: continue
            for ax,(a,b) in enumerate([(r.x0,r.x1),(r.y0,r.y1),(r.z0,r.z1)]):
                if 0.<a<Ls[ax]: g[ax]=insert_layer(g[ax],a*1e-9,0.5e-9,1.3,-1)
                if 0.<b<Ls[ax]: g[ax]=insert_layer(g[ax],b*1e-9,0.5e-9,1.3,+1)
        # interface sheet charges (2DEG): refine on both sides
        for sh in t.get("sheets",[]):
            ax={"x":0,"y":1,"z":2}[sh.axis]
            g[ax]=insert_layer(g[ax],sh.pos*1e-9,0.5e-9,1.3,0)
        # Pin every insulator-region and buried-electrode face with a node on
        # each side (0.5 nm apart).  The half-open ownership rule gives the node
        # ON a plane to the region above it, so without a node just below the
        # far face that face snaps back to the previous grid line: on a coarse
        # grid a 100 nm trench oxide became ~165 nm on one sidewall and 100 nm
        # on the other, and a symmetric IGBT cell conducted 70 % more on one side.
        pinned=[set(),set(),set()]
        def _pin(ax,pos_m):
            arr=np.asarray(g[ax],dtype=float)
            for q in (pos_m-0.5e-9,pos_m):
                if not (arr[0]<q<arr[-1]): continue
                # a pinned node replaces generated nodes closer than 0.05 nm
                # (a graded node 0.01 nm from a pin only degrades conditioning)
                near=(np.abs(arr-q)<0.05e-9)&(np.abs(arr-q)>1e-13)
                near&=~np.isin(arr,list(pinned[ax]))&(arr>arr[0])&(arr<arr[-1])
                arr=np.append(arr[~near],q); pinned[ax].add(q)
            out=np.array(sorted(set(arr.tolist()))); g[ax]=out[np.concatenate(([True],np.diff(out)>1e-13))]
        for r in t.get("regs",[]):
            if r.mat not in INSULATORS: continue
            for ax,(a,b) in enumerate([(r.x0,r.x1),(r.y0,r.y1),(r.z0,r.z1)]):
                for v in (a,b):
                    if 0.<v<Ls[ax]: _pin(ax,v*1e-9)
        for c in t.get("cons",[]):
            if c.face!=6: continue
            for ax,(fa,fb) in enumerate([(c.i0p,c.i1p),(c.j0p,c.j1p),(c.k0p,c.k1p)]):
                for f in (fa,fb):
                    v=f*Ls[ax]
                    if 0.<v<Ls[ax]: _pin(ax,v*1e-9)
        # template-requested boundary layers: (axis, position nm, side)
        for axn,pos,side in (t.get("refine") or []):
            ax={"x":0,"y":1,"z":2}[axn]
            g[ax]=insert_layer(g[ax],pos*1e-9,0.5e-9,1.3,side)
        # no generated node within 0.05 nm of a pinned one
        for ax in range(3):
            if not pinned[ax]: continue
            arr=np.asarray(g[ax],dtype=float); pv=np.array(sorted(pinned[ax]))
            d=np.min(np.abs(arr[:,None]-pv[None,:]),axis=1)
            keep=(d>=0.05e-9)|(d<1e-13)|(arr==arr[0])|(arr==arr[-1])
            g[ax]=arr[keep]
        xs,ys,zs=g
    # Explicit node positions (nm): geometry-critical nodes such as both faces
    # of a 1 nm gate oxide.  Applied on every mesh type; any other node closer
    # than 0.04 nm to one of them is dropped (no near-zero spacings), and so is
    # every generated node inside a "clear" range (a thin gate stack must
    # contain exactly the template's nodes, or its layer thicknesses move).
    if t.get("clear"):
        g=[np.asarray(xs,dtype=float),np.asarray(ys,dtype=float),np.asarray(zs,dtype=float)]
        for axn,rngs in t["clear"].items():
            ax={"x":0,"y":1,"z":2}[axn]; arr=g[ax]
            for a,b in rngs:
                arr=arr[~((arr>a*1e-9+1e-13)&(arr<b*1e-9-1e-13))]
            g[ax]=arr
        xs,ys,zs=g
    if t.get("nodes"):
        g=[np.asarray(xs,dtype=float),np.asarray(ys,dtype=float),np.asarray(zs,dtype=float)]
        for axn,vals in t["nodes"].items():
            ax={"x":0,"y":1,"z":2}[axn]; arr=g[ax]
            keep=[q*1e-9 for q in vals if arr[0]<q*1e-9<arr[-1]]
            if not keep: continue
            kp=np.array(sorted(keep))
            other=[v for v in arr.tolist() if np.min(np.abs(kp-v))>0.04e-9 or v in (arr[0],arr[-1])]
            out=np.array(sorted(set(other)|set(keep)))
            g[ax]=out[np.concatenate(([True],np.diff(out)>1e-13))]
        xs,ys,zs=g
    return xs,ys,zs

def prism_axis(t, axis=None):
    """Prism (extrusion) axis of the triangular mesh: the explicit choice, the
    template's 'extrude' key, the axis of its shapes, else the axis with the
    fewest base nodes (z for a 2-D cross-section)."""
    if axis in ("x","y","z"): return axis
    if t.get("extrude") in ("x","y","z"): return t["extrude"]
    for obj in list(t["regs"])+list(t["cons"]):
        sh=getattr(obj,"shape",None)
        if sh is not None: return sh.axis
    return "xyz"[int(np.argmin([t["Nx"],t["Ny"],t["Nz"]]))]

def owner_grid(regs, xs, ys, zs):
    """Index of the region owning every node of a tensor grid (xs,ys,zs in
    metres), shape (Nz,Ny,Nx): the last region containing it (half-open
    boxes, span_idx, as the core's region fill) cut to its shape; -1 = none."""
    own=np.full((len(zs),len(ys),len(xs)),-1,dtype=np.int32)
    P=None
    for idx,r in enumerate(regs):
        i0,i1=span_idx(xs,r.x0*1e-9,r.x1*1e-9); j0,j1=span_idx(ys,r.y0*1e-9,r.y1*1e-9)
        l0,l1=span_idx(zs,r.z0*1e-9,r.z1*1e-9)
        sub=np.zeros(own.shape,dtype=bool); sub[l0:l1+1,j0:j1+1,i0:i1+1]=True
        if getattr(r,"shape",None) is not None:
            if P is None: P=np.meshgrid(np.asarray(zs)*1e9,np.asarray(ys)*1e9,np.asarray(xs)*1e9,indexing="ij")
            sub&=SM._shape_mask(r.shape,P[2],P[1],P[0],1e-6)
        own[sub]=idx
    return own

def _node_arrays(regs, own, mids):
    """Per-node material index, Nd, Na (m^-3) from the owner indices."""
    own=np.asarray(own).ravel()
    m0=mids[MATS[regs[0].mat]] if regs else 0
    mat=np.full(own.size,m0,dtype=np.int32); Nd=np.zeros(own.size); Na=np.zeros(own.size)
    for idx,r in enumerate(regs):
        sel=own==idx
        if not sel.any(): continue
        mat[sel]=mids[MATS[r.mat]]
        if r.doping>1e11 and r.mat not in INSULATORS:
            if r.dtype=="n": Nd[sel]=r.doping*1e6
            elif r.dtype=="p": Na[sel]=r.doping*1e6
    return mat,Nd,Na

def _ip(a): return np.ascontiguousarray(a,dtype=np.int32).ctypes.data_as(IP)
def _dp(a): return np.ascontiguousarray(a,dtype=np.float64).ctypes.data_as(DP)

def _gate_stack(cid, con, mids, T):
    """Layer stack of a Robin gate for tunnelling (ignored unless tunnelling is on)."""
    lay=con.layers(); ids=[]
    for m,_ in lay:
        mi=MATS.get(m,8)
        if mi not in mids: mids[mi]=api_addmat(mi,T)
        ids.append(mids[mi])
    if ids and min(ids)>=0:
        api_gstack(cid,len(lay),(I*len(lay))(*ids),(D*len(lay))(*[tt*1e-9 for _,tt in lay]))

_PRISM=None          # the prism mesh of the device in the core (None: tensor mesh)

def setup_device(t, solver=(2000,1000,500,1e-8,1e-6,1.,1.), models=7, T=300., graded=True,
                 volt=None):
    """Load a device (template dict) into the C core. volt: {label: V} overrides.
    t["mesh"] = "tri" builds the triangular-prism mesh (axis t["prism"] or
    automatic), anything else the tensor mesh.  Returns the tensor lines
    (xs,ys,zs) in metres - on a prism mesh, the grid results are shown on."""
    global _PRISM
    if t.get("mesh")=="tri":
        return _setup_tri(t,solver,models,T,volt)
    _PRISM=None
    api_init()
    mids={}
    for r in t["regs"]:
        m=MATS[r.mat]
        if m not in mids: mids[m]=api_addmat(m,T)
    xs,ys,zs=make_grid(t,graded)
    Nx,Ny,Nz=len(xs),len(ys),len(zs)
    api_solver(*solver)
    ok=api_gridc(Nx,Ny,Nz,(D*Nx)(*xs),(D*Ny)(*ys),(D*Nz)(*zs))
    if not ok:
        N_total=Nx*Ny*Nz
        raise MemoryError(
            f"Could not allocate {Nx}x{Ny}x{Nz} = {N_total:,} nodes "
            f"(~{N_total*BYTES_PER_NODE/1024**3:.1f} GB).\nTry a coarser grid.")
    api_models(models)
    if any(getattr(r,"shape",None) is not None for r in t["regs"]):
        # shaped regions: node by node (a staircase of the outline)
        mat,Nd,Na=_node_arrays(t["regs"],owner_grid(t["regs"],xs,ys,zs),mids)
        api_ndata(_ip(mat),_dp(Nd),_dp(Na),_dp(np.zeros(len(mat))))
    else:
        if t["regs"]:
            m0=MATS[t["regs"][0].mat]; api_region(0,0,0,Nx-1,Ny-1,Nz-1,mids[m0],0.,0.)
        for r in t["regs"]:
            mx=mids[MATS[r.mat]]
            Nd=r.doping*1e6 if r.dtype=="n" else 0.
            Na=r.doping*1e6 if r.dtype=="p" else 0.
            if r.doping<=1e11 or r.mat in INSULATORS: Nd=Na=0.   # intrinsic / oxide
            # Half-open physical mapping, identical rule to contacts (span_idx)
            ix0,ix1=span_idx(xs,r.x0*1e-9,r.x1*1e-9)
            iy0,iy1=span_idx(ys,r.y0*1e-9,r.y1*1e-9)
            iz0,iz1=span_idx(zs,r.z0*1e-9,r.z1*1e-9)
            api_region(ix0,iy0,iz0,ix1,iy1,iz1,mx,Nd,Na)
    coords={"x":xs,"y":ys,"z":zs}; other={"x":("y","z"),"y":("x","z"),"z":("x","y")}
    for sh in t.get("sheets",[]):
        c=np.asarray(coords[sh.axis]); idx=int(np.argmin(np.abs(c-sh.pos*1e-9)))
        oa,ob=other[sh.axis]
        a0,a1=span_idx(coords[oa],sh.a0*1e-9,sh.a1*1e-9)
        b0,b1=span_idx(coords[ob],sh.b0*1e-9,sh.b1*1e-9)
        api_sheet({"x":0,"y":1,"z":2}[sh.axis],idx,a0,a1,b0,b1,sh.sigma*1e4)
    api_clrcon()
    P3=None
    for con in t["cons"]:
        i0,i1,j0,j1,k0,k1=con.to_idx6(xs,ys,zs)
        V=volt.get(con.label,con.V) if volt else con.V
        robin=con.bc==2 and con.face<=5           # oxide (Robin) gate on a face
        tox=con.eot() if robin else con.tox
        cid=api_addconx(con.face,i0,i1,j0,j1,k0,k1,con.bc,V,con.pm,tox*1e-9,con.nox*1e4)
        if robin and cid>=0: _gate_stack(cid,con,mids,T)
        if cid>=0 and getattr(con,"shape",None) is not None:
            # shaped box contact: its nodes as a list (staircase of the outline)
            if P3 is None: P3=np.meshgrid(np.asarray(zs)*1e9,np.asarray(ys)*1e9,np.asarray(xs)*1e9,indexing="ij")
            L3=(t["Lx"],t["Ly"],t["Lz"])
            m=SM.box_contact_mask(con,P3[2],P3[1],P3[0],L3).ravel()
            nodes=np.nonzero(m)[0]
            if len(nodes): api_cnodes(cid,len(nodes),_ip(nodes),_dp(np.zeros(len(nodes))))
    return xs,ys,zs

def _setup_tri(t, solver, models, T, volt):
    """setup_device on the triangular-prism mesh."""
    global _PRISM
    _PRISM=None
    api_init()
    mids={}
    for r in t["regs"]:
        m=MATS[r.mat]
        if m not in mids: mids[m]=api_addmat(m,T)
    xs,ys,zs=make_grid(t,True)
    pm=SM.PrismMesh(t,(xs,ys,zs),prism_axis(t,t.get("prism")),INSULATORS,thin=t.get("tri_thin",True))
    api_solver(*solver)
    N,E,rp,cj,h,ar,vol,pos=pm.graph()
    if not api_graph(N,E,_ip(rp),_ip(cj),_dp(h),_dp(ar),_dp(vol),_dp(pos)):
        raise MemoryError(f"Could not set up the prism mesh ({N:,} nodes, {E:,} edges).")
    api_models(models)
    own=pm.node_regions()
    mat,Nd,Na=_node_arrays(t["regs"],own,mids)
    api_ndata(_ip(mat),_dp(Nd),_dp(Na),_dp(pm.sheet_charge()))
    ins_id={mids[MATS[m]] for m in INSULATORS if MATS[m] in mids}
    kind=np.where(np.isin(mat,list(ins_id)),1,0).astype(np.int8)
    api_clrcon(); robin=[]
    for con in t["cons"]:
        V=volt.get(con.label,con.V) if volt else con.V
        rob=con.bc==2 and con.face<=5
        tox=con.eot() if rob else con.tox
        cid=api_addconx(con.face,0,0,0,0,0,0,con.bc,V,con.pm,tox*1e-9,con.nox*1e4)
        if cid<0: continue
        nodes,areas=pm.contact_nodes(con)
        api_cnodes(cid,len(nodes),_ip(nodes),_dp(areas))
        if rob:
            _gate_stack(cid,con,mids,T)
            sn=nodes[kind[nodes]==0]
            if len(sn): robin.append((sn,con.face//2,con.face%2))
            kind[nodes[kind[nodes]!=0]]=2
        else:
            kind[nodes[(kind[nodes]!=0)|(con.bc==2)]]=2
    api_walls(_dp(pm.walls(kind,robin)))
    pm.kind_setup=kind; pm.mat_setup=mat
    _PRISM=pm
    return xs,ys,zs

# ══════════════════════════════════════════════════════════════
#  PULL SOLUTION
# ══════════════════════════════════════════════════════════════
_PULL_Q=[("phi",0,1.),("n",1,1e-6),("p",2,1e-6),("Ec",3,1.),("Ev",4,1.),
         ("Efn",5,1.),("Efp",6,1.),("Jx",7,1e-4),("Jy",8,1e-4),("Jz",9,1e-4),
         ("mn",10,1e4),("mp",11,1e4),("R",12,1e-6),("kind",13,1.),("Nnet",14,1e-6),
         ("Jg",27,1e-4),("Lqn",28,1.),("Lqp",29,1.)]  # J in A/cm2, quantum potentials in eV
def pull():
    if _PRISM is not None and api_isgraph():
        d=_pull_prism(_PRISM)
    else:
        Nx,Ny,Nz=api_Nx(),api_Ny(),api_Nz(); N=Nx*Ny*Nz
        if N<=0: return None
        sh=(Nz,Ny,Nx)
        x=np.array([api_xs(i)*1e9 for i in range(Nx)])
        y=np.array([api_ys(j)*1e9 for j in range(Ny)])
        z=np.array([api_zs(l)*1e9 for l in range(Nz)])
        buf=(D*N)()                        # bulk memcpy: one C call per quantity
        d=dict(Nx=Nx,Ny=Ny,Nz=Nz,x=x,y=y,z=z,mesh="rect",nodes=N)
        for key,idx,sc in _PULL_Q:
            api_bulk(idx,buf)
            d[key]=np.frombuffer(buf,dtype=np.float64,count=N).copy().reshape(sh)*sc
    # oxide / metal nodes carry no carriers: blank them (NaN) in carrier plots
    semi = d["kind"]==0
    for key in ("n","p","mn","mp","R","Efn","Efp"):
        d[key]=np.where(semi,d[key],np.nan)
    for key in ("n","p"):
        d[key]=np.where(semi,np.maximum(d[key],1e-30),np.nan)
    d["semi"]=semi
    # terminal currents
    nc=_lib.api3_get_ncon() if hasattr(_lib,"api3_get_ncon") else 0
    d["I"]=[api_conI(c) for c in range(nc)]; d["A"]=[api_conA(c) for c in range(nc)]
    d["Ifl"]=[api_conF(c) for c in range(nc)]
    d["Vs"]=[api_Vsol(c) for c in range(nc)]
    return d

def _pull_prism(pm):
    """Solution of a prism mesh resampled onto its tensor lines (every node of
    the prism mesh is one of them; the thinned-out ones are interpolated
    inside their triangle - from semiconductor nodes only for carrier
    quantities, from the dominant node's material for band edges)."""
    N=pm.N; buf=(D*N)()
    def get(idx,sc):
        api_bulk(idx,buf); return np.frombuffer(buf,dtype=np.float64,count=N).copy()*sc
    kind=get(13,1.); semi=kind==0
    mat=getattr(pm,"mat_setup",None)
    x,y,z=(c*1.0 for c in pm.lines)
    d=dict(Nx=len(x),Ny=len(y),Nz=len(z),x=x,y=y,z=z,mesh="tri",nodes=N,prism=pm)
    raw={"kind":kind,"semi":semi}
    for key,idx,sc in _PULL_Q:
        if key=="kind": d[key]=pm.to_display(kind,"dom"); continue
        v=get(idx,sc); raw[key]=v
        if key in ("n","p"): d[key]=pm.to_display(np.maximum(v,1e-300),"log",mask=semi)
        elif key in ("Nnet",): d[key]=pm.to_display(v,"dom")
        elif key in ("Ec","Ev"): d[key]=pm.to_display(v,"lin",mat=mat)
        elif key=="phi": d[key]=pm.to_display(v,"lin")
        else: d[key]=pm.to_display(v,"lin",mask=semi)
    # the node values themselves (cross-section plots on the triangulation)
    for key in ("n","p","mn","mp","R","Efn","Efp"): raw[key]=np.where(semi,raw[key],np.nan)
    for key in ("n","p"): raw[key]=np.where(semi,np.maximum(raw[key],1e-30),np.nan)
    d["raw"]=raw
    return d
_lib.api3_get_ncon.restype=I; _lib.api3_get_ncon.argtypes=[]

# ══════════════════════════════════════════════════════════════
#  EXPERIMENT ENGINE  (GUI-independent: nets, sweeps, floating nodes)
# ══════════════════════════════════════════════════════════════
def _opt(fn,res,*args,default=None):
    """Optional binding: a no-op returning `default` if the .so predates it."""
    try: return _f(fn,res,*args)
    except AttributeError: return (lambda *a: default)
api_stop    = _opt("api3_request_stop",V)          # cooperative stop (any thread)
api_unstop  = _opt("api3_clear_stop",V)
api_prog    = _opt("api3_get_progress",V,DP)       # 8 doubles, see core
api_Vsol    = _opt("api3_get_contact_Vsolved",D,I,default=float("nan"))
api_threads = _opt("api3_set_threads",V,I)
api_aa      = _opt("api3_get_aa_depth",I,default=0)
api_gstack  = _opt("api3_set_gate_stack",I,I,I,IP,DP,default=0)   # cid,nl,mat idx,t(m)
api_ntp     = _opt("api3_get_ntp",I,I,default=0)                   # tunnelling paths (-1: all)
# general box-method meshes (triangular prisms) and node-wise device data
api_graph   = _opt("api3_set_mesh_graph",I,I,I,IP,IP,DP,DP,DP,DP,default=0)  # N,E,rp,cj,h,ar,vol,pos
api_isgraph = _opt("api3_is_graph",I,default=0)
api_ndata   = _opt("api3_set_node_data",V,IP,DP,DP,DP)                     # mat,Nd,Na,Nf per node
api_cnodes  = _opt("api3_set_con_nodes",I,I,I,IP,DP,default=0)            # cid,n,nodes,areas
api_walls   = _opt("api3_set_walls",V,DP)                                  # 6 per node
# model switches (api3_set_models): bit 0..5
M_BGN,M_TAU,M_VSAT,M_QC,M_LB,M_TUN = 1,2,4,8,16,32

def progress():
    """(t of the running ramp, Gummel it, residual, steps, phase, total its)"""
    b=(D*8)(); api_prog(b); return tuple(b[:6])

def con_groups(cons):
    """Electrical nodes of the device: [(key, [contact indices])].  Contacts
    that share a net name form one node (key = net name); every other contact
    is its own node (key = its label; made unique if labels repeat)."""
    out=[]; where={}
    for k,c in enumerate(cons):
        net=(getattr(c,"net","") or "").strip()
        if net:
            if net in where: out[where[net]][1].append(k); continue
            where[net]=len(out); out.append([net,[k]])
        else:
            key=c.label; n=2
            while key in where: key=f"{c.label} #{n}"; n+=1
            where[key]=len(out); out.append([key,[k]])
    return [(k,v) for k,v in out]

def group_of(cons,key):
    """Contact indices of the node `key` (net name or label), [] if unknown."""
    for k,idx in con_groups(cons):
        if k==key: return list(idx)
    return []

def group_label(cons,key,idx):
    if len(idx)>1 or getattr(cons[idx[0]],"net",""):
        return f"{key}  [net: {len(idx)} contact{'s' if len(idx)>1 else ''}]"
    return key

class Job:
    """Progress/stop record shared between a worker thread and the GUI."""
    def __init__(self,kind="run",N=1):
        self.kind=kind; self.N=max(1,N); self.k=0; self.ev=0; self.stop=False
        self.msg=""; self.t0=time.time()
    def request_stop(self):
        self.stop=True; api_stop()

def _fl_eval(fidx,V):
    """Put the floating node at V and solve: returns (I_net, tol, ok, V_reached).
    I_net = total current INTO the device through the node's contacts."""
    for c in fidx: api_setV(c,V)
    ok=api_resolve()
    Vr=api_Vsol(fidx[0])
    if Vr!=Vr: Vr=V
    nc=_lib.api3_get_ncon()
    Inet=sum(api_conI(c) for c in fidx)
    Imax=max([abs(api_conI(c)) for c in range(nc)]+[0.])
    Fl=sum(api_conF(c) for c in fidx)
    return Inet, max(1e-3*Imax, 3.*Fl, 1e-30), bool(ok), Vr

def solve_floating(fidx,V0,lo,hi,job=None,slope=None,vtol=1e-9,maxev=40):
    """Zero-current (floating) node: find V in [lo,hi] such that the total
    current into the device through the contacts `fidx` vanishes.
    The current is non-decreasing in V and, by the maximum principle, the root
    lies between the lowest and highest applied voltages (lo, hi), so the
    bracket needs no evaluation of its ends.  Safeguarded secant (Illinois
    regula falsi once both signs are known, bisection as a fallback); every
    evaluation continues from the previous solution.
    Returns (V, I_net, converged, evaluations, dI/dV estimate)."""
    a,b=float(lo),float(hi); fa=fb=None; side=0
    x=min(max(float(V0),a),b); hist=[]; ok=False; fx=0.; xr=x; found=False
    for ev in range(1,maxev+1):
        if job is not None:
            if job.stop: break
            job.ev=ev
        fx,tol,ok,xr=_fl_eval(fidx,x)
        hist.append((xr,fx))
        if ok and abs(fx)<=tol: found=True; break
        if fx>0.: b,fb=xr,fx; side=(side-1 if side<0 else -1)
        else:     a,fa=xr,fx; side=(side+1 if side>0 else 1)
        if b-a<=vtol: found=ok; break
        cand=None
        if len(hist)>=2:                            # secant through the last two points
            (x1,f1),(x2,f2)=hist[-2],hist[-1]
            if f2!=f1 and x2!=x1: cand=x2-f2*(x2-x1)/(f2-f1)
        elif slope and slope>0.: cand=xr-fx/slope   # slope from the previous bias point
        if fa is not None and fb is not None and not (cand is not None and a<cand<b):
            fa2,fb2=fa,fb                           # Illinois: damp the stale end
            if side>=2: fb2=0.5*fb                  # a replaced twice -> b is stale
            if side<=-2: fa2=0.5*fa                 # b replaced twice -> a is stale
            cand=b-fb2*(b-a)/(fb2-fa2) if fb2!=fa2 else 0.5*(a+b)
        if cand is None or not (a<cand<b):
            cand=0.5*(a+b)
        # keep a minimum step so the secant cannot stall on a flat plateau
        if abs(cand-xr)<0.25*vtol: cand=xr+(0.25*vtol if fx<0 else -0.25*vtol)
        # The node sits on a rail when one transistor is off (inverter output
        # within nV of VDD): the current there is a steep linear function of
        # V, so the tolerance is on the CURRENT - a bracket of 0.1 mV would
        # leave microamps of spurious rail-to-rail current.
        x=min(max(cand,a),b)
    sl=None
    if len(hist)>=2:
        (x1,f1),(x2,f2)=hist[-2],hist[-1]
        if x2!=x1 and (f2-f1)/(x2-x1)>0: sl=(f2-f1)/(x2-x1)
    return xr,fx,found,len(hist),(sl or slope)

def run_sweep(cons,tgt,Vs,Ve,Np,flt=None,job=None,on_point=None,vapp=None,refine=True):
    """Sweep the node `tgt` (contact indices, driven together) from Vs to Ve in
    Np points, continuing from the present C state (setup_device first).
    vapp: applied voltage of every contact (default: the Con.V values); the
    non-swept, non-floating contacts are held there.
    flt: contact indices of a floating node solved for zero current at every
    point (e.g. an inverter output).  With refine=True an output jump larger
    than 8 % of the voltage span between neighbouring points is bisected
    (the steep part of a transfer curve gets its own points).
    on_point(k, rec) is called after every solved point (calling thread).
    Returns a dict of numpy arrays ordered along the sweep."""
    nc=len(cons); Np=max(1,int(Np))
    vapp=list(vapp) if vapp is not None else [c.V for c in cons]
    tgt=list(tgt); flt=[c for c in (flt or []) if c not in tgt]
    for c in range(nc):
        if c not in tgt and c not in flt: api_setV(c,vapp[c])
    Vv=np.linspace(Vs,Ve,Np) if Np>1 else np.array([float(Vs)])
    R=dict(V=[],I=[],F=[],conv=[],Vf=[],nev=[],t=[])
    fixed=[vapp[c] for c in range(nc) if c not in flt and c not in tgt]
    span=max(fixed+[Vs,Ve])-min(fixed+[Vs,Ve]) or 1.
    st={"slope":None}
    if job is not None: job.N=Np; job.extra=0
    def solve_at(Vt,guess):
        for c in tgt: api_setV(c,Vt)
        if flt:
            f=fixed+[Vt]; lo,hi=min(f),max(f)
            Vx,fx,ok,ne,sl=solve_floating(flt,min(max(guess,lo),hi),lo,hi,job,st["slope"])
            st["slope"]=sl
        else:
            ok=bool(api_resolve()); ne=1; Vx=float("nan")
        Vr=api_Vsol(tgt[0])
        rec=dict(V=(Vr if Vr==Vr else Vt),I=[api_conI(c) for c in range(nc)],
                 F=[api_conF(c) for c in range(nc)],conv=bool(ok),Vf=Vx,nev=ne,t=time.time())
        for kk in R: R[kk].append(rec[kk])
        if on_point: on_point(len(R["V"])-1,rec)
        return rec
    def guess_at(Vt):
        if not flt: return 0.
        if not R["V"]: return vapp[flt[0]]
        V=np.array(R["V"]); Vf=np.array(R["Vf"],dtype=float)
        o=np.argsort(np.abs(V-Vt))
        if len(o)>=2 and V[o[0]]!=V[o[1]]:
            a,b=o[0],o[1]; return Vf[a]+(Vf[a]-Vf[b])*(Vt-V[a])/(V[a]-V[b])
        return Vf[o[0]]
    dVmax=0.08*span; dVin_min=abs(Ve-Vs)/max(Np-1,1)/32.; budget=3*Np
    prev=None
    for k,Vt in enumerate(Vv):
        if job is not None:
            if job.stop: break
            job.k=k; job.ev=0
        rec=solve_at(Vt,guess_at(Vt))
        if flt and refine and prev is not None:
            stack=[(prev,rec)]
            while stack and budget>0 and not (job is not None and job.stop):
                a,b=stack.pop()
                if abs(b["Vf"]-a["Vf"])<=dVmax or abs(b["V"]-a["V"])<=dVin_min: continue
                Vm=0.5*(a["V"]+b["V"])
                m=solve_at(Vm,0.5*(a["Vf"]+b["Vf"])); budget-=1
                if job is not None: job.extra=getattr(job,"extra",0)+1
                stack.append((m,b)); stack.append((a,m))
        prev=rec
    out={k:np.array(v) for k,v in R.items()}
    if len(out["V"]):
        o=np.argsort(out["V"],kind="stable")
        if Ve<Vs: o=o[::-1]
        out={k:v[o] for k,v in out.items()}
    out["conv"]=out["conv"].astype(bool) if len(out["conv"]) else np.zeros(0,bool)
    return out

def run_point(cons,flt=None,job=None,vapp=None):
    """Solve at the applied voltages (continuing from the C state); with a
    floating node, solve it for zero current.  Returns (ok, V_float)."""
    nc=len(cons); vapp=list(vapp) if vapp is not None else [c.V for c in cons]
    flt=list(flt or [])
    for c in range(nc):
        if c not in flt: api_setV(c,vapp[c])
    if not flt: return bool(api_resolve()), None
    fixed=[vapp[c] for c in range(nc) if c not in flt]
    lo,hi=(min(fixed),max(fixed)) if fixed else (-1.,1.)
    Vx,fx,ok,ne,_=solve_floating(flt,vapp[flt[0]],lo,hi,job)
    return ok,Vx

def _shape_key(o):
    sh=getattr(o,"shape",None)
    return None if sh is None else repr(sh.to_dict())

def _mesh_desc(s):
    """'Nx×Ny×Nz nodes' or the prism mesh size of a pulled solution."""
    if s.get("mesh")=="tri":
        pm=s.get("prism")
        return (f"{s['nodes']:,} nodes (prisms {pm.N2:,} × {pm.NL})" if pm is not None else f"{s['nodes']:,} nodes")
    return f"{s['Nx']}×{s['Ny']}×{s['Nz']} nodes"

def device_signature(t,solver,models,T,graded):
    """Everything that requires a fresh setup_device (not the voltages)."""
    regs=[(r.mat,r.x0,r.y0,r.z0,r.x1,r.y1,r.z1,r.dtype,r.doping,_shape_key(r)) for r in t["regs"]]
    cons=[(c.face,c.i0p,c.i1p,c.j0p,c.j1p,c.k0p,c.k1p,c.bc,c.pm,c.tox,c.nox,repr(c.stack),_shape_key(c)) for c in t["cons"]]
    shs=[(s.axis,s.pos,s.a0,s.a1,s.b0,s.b1,s.sigma) for s in t.get("sheets",[])]
    return repr((t["Lx"],t["Ly"],t["Lz"],t["Nx"],t["Ny"],t["Nz"],regs,cons,shs,
                 sorted(t.get("dense_x",[])),sorted(t.get("dense_y",[])),sorted(t.get("dense_z",[])),
                 repr(t.get("nodes")),repr(t.get("clear")),repr(t.get("refine")),tuple(solver),models,T,graded,
                 t.get("mesh","rect"),t.get("prism"),t.get("extrude")))

# ══════════════════════════════════════════════════════════════
#  GUI TOOLKIT: fonts, two-way scrolling pages, tooltips, help text
# ══════════════════════════════════════════════════════════════
import queue, json, copy, textwrap
import tkinter.font as tkfont
from matplotlib.ticker import FuncFormatter

VERSION="7.4"
SETTINGS_FILE=os.path.join(os.path.expanduser("~"),".semisim3d.json")

def _first_family(root,cands,fallback):
    try: fams=set(tkfont.families(root))
    except Exception: fams=set()
    for c in cands:
        if c in fams: return c
    return fallback

class ScrollFrame(tk.Frame):
    """A page whose content (.inner) scrolls vertically AND horizontally.
    The inner frame is stretched to the viewport when it is smaller, so
    widgets packed with expand=True still fill the page."""
    def __init__(self,parent,bg,**kw):
        super().__init__(parent,bg=bg,**kw)
        self.canvas=tk.Canvas(self,bg=bg,highlightthickness=0,bd=0,
                              xscrollincrement=24,yscrollincrement=24)
        self.vsb=ttk.Scrollbar(self,orient="vertical",command=self.canvas.yview)
        self.hsb=ttk.Scrollbar(self,orient="horizontal",command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=self.vsb.set,xscrollcommand=self.hsb.set)
        self.canvas.grid(row=0,column=0,sticky="nsew")
        self.vsb.grid(row=0,column=1,sticky="ns"); self.hsb.grid(row=1,column=0,sticky="ew")
        self.rowconfigure(0,weight=1); self.columnconfigure(0,weight=1)
        self.inner=tk.Frame(self.canvas,bg=bg)
        self._win=self.canvas.create_window(0,0,window=self.inner,anchor="nw")
        self.inner.bind("<Configure>",self._sync); self.canvas.bind("<Configure>",self._sync)
        self.canvas.wheel_target=True
    def _sync(self,_=None):
        rw,rh=self.inner.winfo_reqwidth(),self.inner.winfo_reqheight()
        cw,ch=self.canvas.winfo_width(),self.canvas.winfo_height()
        W,H=max(rw,cw),max(rh,ch)
        self.canvas.itemconfigure(self._win,width=W,height=H)
        self.canvas.configure(scrollregion=(0,0,W,H))

class Tip:
    """Hover tooltip."""
    def __init__(self,w,text,font=None):
        self.w=w; self.text=text; self.tw=None; self.id=None; self.font=font
        w.bind("<Enter>",self._sched,add="+"); w.bind("<Leave>",self._hide,add="+")
        w.bind("<ButtonPress>",self._hide,add="+")
    def _sched(self,_=None):
        self._cancel(); self.id=self.w.after(650,self._show)
    def _cancel(self):
        if self.id:
            try: self.w.after_cancel(self.id)
            except Exception: pass
            self.id=None
    def _show(self):
        self.id=None
        try:
            x=self.w.winfo_rootx()+12; y=self.w.winfo_rooty()+self.w.winfo_height()+4
            self.tw=tk.Toplevel(self.w); self.tw.wm_overrideredirect(True)
            self.tw.wm_geometry(f"+{x}+{y}")
            tk.Label(self.tw,text=self.text,bg="#16263b",fg="#e6f0ff",relief="solid",bd=1,
                     justify="left",padx=7,pady=4,wraplength=420,font=self.font).pack()
        except Exception: self.tw=None
    def _hide(self,_=None):
        self._cancel()
        if self.tw is not None:
            try: self.tw.destroy()
            except Exception: pass
            self.tw=None

def render_markup(txt,body):
    """Tiny markup -> Tk text tags:  '# ' heading, '## ' subheading,
    '- ' bullet, '  - ' sub-bullet, '  ' (two spaces) preformatted line,
    `code` spans."""
    import re
    for line in body.strip("\n").split("\n"):
        if line.startswith("## "):
            txt.insert("end",line[3:]+"\n","h2"); continue
        if line.startswith("# "):
            txt.insert("end",line[2:]+"\n","h1"); continue
        tag=()
        if line.startswith("- "):
            txt.insert("end","  •  ","bullet"); line=line[2:]; tag=("bullet",)
        elif line.startswith("  - "):
            txt.insert("end","      ◦  ","bullet2"); line=line[4:]; tag=("bullet2",)
        elif line.startswith("  "):
            txt.insert("end",line+"\n","pre"); continue
        parts=re.split(r"(`[^`]+`)",line)
        for p in parts:
            if p.startswith("`") and p.endswith("`") and len(p)>1:
                txt.insert("end",p[1:-1],("code",)+tag)
            else:
                txt.insert("end",p,tag)
        txt.insert("end","\n",tag)

def help_topics():
    """[(title, markup)] - the user manual shown in Help (F1)."""
    cats="\n".join(
        f"## {cat}\n"+"\n".join(f"- `{n}` - {TEMPLATES[n]().get('desc','')}" for n in names)
        for cat,names in template_groups())
    return [
("Quick start", """
# Quick start
- Pick a device from the `Device ▾` menu in the toolbar (or `Devices` in the menu bar). The structure appears in the `Structure` tab.
- Press `▶ Run` (F5) to solve the device at the voltages listed in the editor's `Contacts` page.
- Press `↔ Sweep` (F6) to run the I-V sweep defined in the `Sweep` page. The curve builds up live in the `I-V` tab.
- Press `■ Stop` (Esc) at any time. The solver keeps the last converged bias point, so all plots stay valid.
- Open or close the editor panel with `☰` (Ctrl+E). With the panel closed the plots use the whole window.
- Every panel scrolls both ways. Use the scroll bars, the mouse wheel, or Shift + wheel for sideways scrolling.
- Text too small or too large? Use View → Text size, or Ctrl + / Ctrl − (Ctrl 0 resets). Plot labels scale too.

Typical workflow for a new design: load the closest template, change regions and contacts, press Run, inspect the Bands/Current/Field tabs, then sweep.
"""),
("Window layout", """
# Window layout
## Toolbar
- `☰` shows or hides the editor panel (drawer). Drag the divider to resize it.
- `Device ▾` opens the template library by category. `Browse…` shows full descriptions.
- `▶ Run`, `↔ Sweep` and `■ Stop` control the solver. The progress bar shows the bias ramp (Run) or the sweep points (Sweep). The text beside it shows the continuation step, the Gummel iteration and the residual.
## Editor panel (drawer) pages
- `Device`: the template description, domain size, base mesh, mesh type (rectangular or triangular prisms) and mesh statistics.
- `Regions`: the material/doping blocks. Later rows override earlier ones where they overlap.
- `Contacts`: electrodes, their voltages, work functions, gate oxides and nets.
- `Solver`: temperature, physical models, iteration limits, tolerances and CPU threads.
- `Sweep`: the swept terminal, the voltage range, a floating (zero-current) terminal and the plotted current.
## Plot tabs
- `Structure`: 3-D view and a cross-section with the real mesh (lines, or the triangles of a prism mesh). Pick the plane XY/YZ/XZ and its position.
- `3D Slicer`: any solution quantity on three orthogonal planes.
- `Bands`: Ec, Ev and the quasi-Fermi levels along a line cut, plus n and p.
- `Current`: current density along a cut, recombination, the in-plane current field and a 3-D arrow field.
- `I-V`: the sweep result. For a transfer curve (floating output) it shows Vout(Vin), the gain and the supply current.
- `Field & Mobility`: electric field and carrier mobilities along a cut.
- `Results`: a table of every terminal (voltage, current, current density, resolution), extracted parameters and a run log.
Each plot has the matplotlib toolbar (zoom, pan, save) under it.
"""),
("Device templates", """
# Device templates
The library has textbook structures and `RW:` devices modelled on real commercial parts. RW devices use physically consistent dopings and thicknesses for their voltage class, not reverse-engineered dies. Currents are per simulated cell (a few µm wide, Lz deep). Real parts contain thousands of cells in parallel.

"""+cats+"""

## CMOS inverter
Two complementary MOSFETs share the gate net `Vin` and the drain net `Out`.
The default sweep is the voltage transfer curve. `Vin` goes from 0 to 2.5 V, and at every point `Out` is solved for zero current (an unloaded output).
The I-V tab reports:
- the switching threshold V_M, where Vout = Vin;
- the maximum gain;
- V_OH, V_OL, V_IL and V_IH, the points where the gain equals −1;
- the noise margins NM_L = V_IL − V_OL and NM_H = V_OH − V_IH;
- the short-circuit supply-current peak.
Both transistors have the same width, so V_M sits below VDD/2 (the PMOS is weaker).
## FinFET
A 10 nm fin wraps the gate on three sides. The three box electrodes form the net `Gate` and are swept as one terminal.
The gate dielectric is a real high-k stack, meshed node by node: 0.5 nm SiO2 interfacial layer + 2.5 nm HfO2 (EOT 0.94 nm).
The structure tab opens on the YZ cross-section through the gate, so the tri-gate and the stack are visible.
The I-V tab reports the subthreshold swing SS, V_th (maximum-gm extrapolation) and I_on/I_off; the semi-log panel also shows the gate tunnelling current.
## Gate dielectrics & tunnelling
- `Gate leakage: 1.2nm SiO2 NMOS` and `Gate leakage: HfO2 high-k NMOS` are the same transistor (EOT 1.2 nm) with two gate stacks. The sweep plots the gate current: about 1e2 A/cm² through 1.2 nm SiO2 at +1 V against about 1e-2 A/cm² through 0.5 nm SiO2 + 3.95 nm HfO2.
- `Fowler-Nordheim: 8nm SiO2 MOS`: gate current through a thick oxide at 4-10 V, where the barrier turns triangular.
"""),
("Regions, mesh and geometry", """
# Regions, mesh and geometry
Coordinates are in nm. x runs left to right. y is depth: y = 0 is the top surface, so the plots show y growing downwards. z is the third direction. Planar devices are thin in z (a 2-D cross-section).
## Regions
A region is a box of one material with one doping type and concentration (cm⁻³). Regions are applied in table order, and a later region overwrites an earlier one where they overlap. So start with a background block and carve wells, sources and drains on top. Use ▲/▼ to reorder.
Insulators (SiO2, Al2O3, HfO2, Si3N4) carry no carriers. They hold fixed potential gradients and shape gate fields; with gate tunnelling on, carriers can cross them into a gate electrode.
## Mesh
Nx, Ny and Nz are the base node counts. With `Graded mesh` on, nodes cluster at every junction and material interface, and the mesher adds:
- boundary layers (0.5 nm first cell) under oxide-gated faces and at insulator interfaces;
- two-sided layers at interface sheet charges (2DEGs);
- pinned nodes on both sides of every oxide face and buried electrode, so thin oxides keep their exact thickness.
The `Device` page shows the final node count and memory estimate (about 760 bytes per node). There is no fixed node limit. The limit is your RAM and your patience: 3-D devices with 10⁵–10⁶ nodes take minutes per bias point.
## Mesh types (Device page)
- `Rectangular`: the tensor mesh above - every x, y and z line runs through the whole device.
- `Triangular (prisms)`: the cross-section normal to the `Prism axis` is a Delaunay triangulation, extruded along the axis over the same 1-D grid the rectangular mesh uses on that axis. The control volumes are the Voronoi cells of the triangles, so the box method (Scharfetter-Gummel, the same physics) runs unchanged. `auto` picks the axis of the template's outlines, else the axis with the fewest base nodes (z for a 2-D cross-section).
  - Box geometry: near every interface, junction, gate and contact edge the triangular mesh uses the rectangular mesh's own nodes, so interfaces and inversion layers are resolved identically; the fine lines that the rectangular mesh drags through the whole device are dropped away from where they are needed (typically 10-50 % fewer nodes, same results within ~1 %).
  - Round and sloped geometry (regions or gates with an outline, see below): the outline is meshed as it is. A node pair straddles every interface (the face between the two nodes IS the interface), a boundary layer grows into the semiconductor, the metal surface of a wrapped gate carries its own row, and the layers of a gate stack line up along the normals, so tunnelling paths run straight through the stack. The rectangular mesh can only draw such a device as a staircase.
## Outlines (shapes)
A region or a Box contact may carry an outline in the plane normal to one axis: a convex polygon with rounded corners, or a circle. The region is then its box cut to the outline (a gate electrode: its box minus the outline). Outlines come from the templates (`FinFET (tapered, rounded fin)`, `GAA nanowire FET (round)`) and are marked ◯ in the tables; editing the box keeps the outline. A gate stack is built as offsets of one outline (constant film thickness, corners rounded r + t), which is what lets the mesher align its layers.
"""),
("Contacts, gates and nets", """
# Contacts, gates and nets
## Faces and ranges
A contact sits on one of the six domain faces or is a `Box` (a buried electrode).
The i/j ranges are percentages of the face's two in-plane axes. The Contacts page shows which axes they are (for example `i→X j→Z` on the Y faces). A Box uses the percentages of X, Y and Z.
## Boundary conditions
- `Ohmic`: charge-neutral, equilibrium carrier densities at the applied voltage. It works over n and p regions at once, like a source contact that also shorts the body tap.
- `Schottky`: metal with work function φm (eV). The barrier is φm − χ for electrons.
- `Gate`: on a semiconductor face, an oxide of thickness tox (EOT, nm) with fixed charge Nox (q/cm²) is modelled as a boundary condition. On insulator nodes or as a Box, it is a metal electrode.
- `Gate stack` (face gates): the layers between the semiconductor and the metal, semiconductor side first, in nm, e.g. `SiO2 0.5, HfO2 2.5`. The gate capacitance uses the stack's EOT (tox is then computed) and gate tunnelling sees every layer. Empty means one SiO2 layer of thickness tox. A Box gate needs no stack: tunnelling follows the meshed insulator regions between the semiconductor and the metal.
Common work functions: n+ poly 4.05–4.1 eV, p+ poly 5.1–5.2 eV, TiN 4.5–4.6 eV, Ni 5.1 eV.
## Nets
Contacts with the same `Net` name are wired together outside the simulated cell. They always carry the same voltage, they are swept together, and their currents add up in the plots. Editing the voltage of one member updates the whole net.
Examples: the NMOS and PMOS gates of an inverter (`Vin`), the three sides of a FinFET gate (`Gate`), and a source and substrate contact tied to ground.
## Sign convention
A terminal current is positive when current flows INTO the device through that terminal. The currents of all terminals add up to zero (Kirchhoff).
"""),
("Running and sweeps", """
# Running and sweeps
## Run (F5)
Solves the device at the contact voltages. If only voltages changed since the last solve, the solver continues from the previous solution, which is much faster. Turn this off in the `Solver` page to always start from equilibrium.
The bias is ramped adaptively from the previous state. MOS gates are ramped first, then the other terminals, like a real measurement. If a step fails, the step shrinks. If it still fails, the last converged bias point is kept and the status says NOT CONVERGED.
## Sweep (F6)
- `Swept terminal`: a contact or a whole net.
- `Start / End / Points`: the voltage range.
- `Floating terminal (I = 0)`: optional. At every point this terminal is solved for zero current, like the output of a logic gate with no load. This gives transfer curves. Steep parts of the curve are refined automatically when `Adaptive refinement` is on.
- `Plot current of`: the terminal whose current appears in the linear panel. `auto` uses the swept terminal, or the terminal with the largest current if the swept one is an insulated gate.
- `BV knee`: a current threshold (A/cm²) for the heuristic breakdown-knee marker. 0 turns it off.
## Stop (Esc)
Stops within one Gummel iteration. A sweep keeps the points already solved, and the device keeps its last converged state.
## Extracted parameters
The I-V tab and the Results tab report what fits the sweep:
- diode ideality factor and J0 (forward exponential branch);
- subthreshold swing, V_th and I_on/I_off (gate sweeps);
- V_M, gain and noise margins (transfer curves);
- a current-knee marker for blocking sweeps.
"""),
("Plots", """
# Plots
- `Structure`: choose the cross-section plane (XY, YZ, XZ) and slide its position. Mesh lines are the real nodes. Hover over the section to read the region, material and doping under the cursor.
- Every cut control has an axis selector and two position sliders for the other two axes. The label shows the position in nm.
- `Bands`: Ec and Ev (solid), the quasi-Fermi levels Efn and Efp (dashed) and a combined Ef, weighted towards the majority carrier. The energy zero is the vacuum level at φ = 0, so absolute values depend on the gauge but differences are physical. The n and p overlay uses a log scale on the right axis.
- `3D Slicer`: pick a quantity and slide the plane. The three small maps are the XY, YZ and XZ planes through the slice point, drawn on the real (non-uniform) mesh.
- `Current`: the in-plane map shows |J| with direction arrows. Arrow length is log-scaled so weak and strong regions are both visible.
- Oxide and metal nodes carry no carriers. They are blank in carrier plots.
- Triangular (prism) mesh: line cuts and maps use the rectangular lines as a display grid (every prism node lies on it or is interpolated inside its triangle, from the same material). In the 3D Slicer the map of the cross-section normal to the prism axis is drawn on the triangulation itself, and the Structure tab's mesh overlay shows the triangles.
- Use the matplotlib toolbar under each plot to zoom and pan. File → Export plot saves the visible tab (PNG/PDF/SVG).
"""),
("Physics and numerics", """
# Physics and numerics
## Equations
Drift-diffusion: Poisson's equation plus electron and hole continuity with Scharfetter-Gummel fluxes, on a non-uniform 3-D box mesh (finite volumes).
## Models
- Heterojunctions (electron affinity rule).
- Fermi statistics are not used (Boltzmann).
- Klaassen band-gap narrowing.
- Arora doping-dependent mobility with velocity saturation.
- SRH with doping-dependent lifetimes, Auger and radiative recombination.
- Oxide gates with fixed charge, interface sheet charges (GaN polarisation) and buried metal electrodes.
## Nanoscale MOS models (Solver page; MOS templates switch them on)
- Quantum confinement, MLDA (modified local-density approximation): at every semiconductor/insulator interface the density of states is cut by the factor 1 − exp(−(z/λ)²), λ = ħ/√(2 m kT), per valley: Si electrons in the Δ2 (m = 0.916, λ = 1.27 nm, 1/3 of the states) and Δ4 valleys (m = 0.19, λ = 2.79 nm, 2/3), holes heavy/light. The inversion layer moves ~1 nm off the oxide: V_th rises by 40-100 mV and the gate capacitance drops. The Bands tab draws Ec + Λn and Ev − Λp (the quantum potentials).
- Lombardi surface mobility (Si): 1/µ = 1/µ_bulk + D/µ_ac + D/µ_sr with the field normal to the current, D = exp(−d/10 nm). Reproduces the universal mobility curve (~400 cm²/Vs at 0.5 MV/cm, ~250 at 1 MV/cm for electrons).
- Gate tunnelling: Tsu-Esaki current with WKB transmission through the gate stack for electrons (conduction band) and holes (valence band), both directions (semiconductor → gate and gate → semiconductor). Direct tunnelling and Fowler-Nordheim are the same formula. Barriers from the electron affinities: SiO2 3.1 eV (tunnelling mass 0.42), HfO2 1.5 eV (0.18), Al2O3 2.7 eV (0.35), Si3N4 2.0 eV (0.5). The supply of an inversion or accumulation layer comes from its pressure on the interface, so the result does not depend on the mesh or on MLDA. The gate current appears in the Results table and the I-V semi-log panel; the 3D Slicer shows log10|J gate| at the interface nodes.
## Solver
- Newton-Poisson with multigrid-preconditioned CG.
- Continuity with BiCGSTAB and aggregation multigrid (on a prism mesh: strength-based pairwise aggregation of the mesh graph, aggregates of up to 8 nodes along the strongest couplings).
- A Gummel map on (ψ, φn, φp) accelerated by Anderson mixing.
- A floating-region balance step for regions tied to no contact.
- Adaptive bias continuation.
- Terminal currents are exact in the discrete Kirchhoff sense and come with a round-off resolution floor. Currents below the floor are shown hollow or omitted.
## Validation
`validate.py` compares the core with closed-form theory: built-in voltage, depletion width, Shockley diode current, SiC blocking field, heterojunction offsets, MOS threshold, MOSFET linear current, the BJT Gummel integral and more. Section 13 checks the triangular-prism mesh: with every rectangular node kept it reproduces the rectangular solution exactly, the thinned mesh agrees within 1 % on MOSFETs, a GaN HEMT and the FinFET, and a round MOS capacitor (Si wire in an oxide shell) matches the radial Poisson-Boltzmann solution within a few mV where a staircase mesh misses by 20-80 mV. Sections 11-12 check the tunnelling kernel against direct integration, gate leakage against textbook values (1.2 nm SiO2: 1e2-1e3 A/cm² at 1 V; ~5 decades per nm), the Fowler-Nordheim slope against measured Si/SiO2 data, the MLDA threshold shift and inversion-layer centroid against quantum-mechanical results, and Lombardi against the universal mobility curve.
"""),
("Limitations", """
# Limitations and troubleshooting
## Not modelled
- impact ionisation: breakdown voltage is NOT predicted, and the knee marker is a heuristic;
- band-to-band tunnelling: Zener behaviour and GIDL are absent;
- optical generation: solar cells show dark I-V only;
- quantisation beyond MLDA: no subband energies, and MLDA underestimates the dark space at very high fields; no quantum confinement in heterostructure 2DEGs;
- gate tunnelling: no image-force barrier lowering, no trap-assisted tunnelling, no valence-electron (EVB) tunnelling;
- thermionic emission at heterointerfaces;
- ballistic transport, strain and remote-phonon scattering in high-k stacks: nanoscale on-currents are drift-diffusion estimates;
- self-heating;
- transients and AC;
- general 3-D (tetrahedral) meshes: the triangular mesh is a 2-D triangulation extruded along one axis, so an outline must be normal to the prism axis (others fall back to a staircase), and geometry along the axis stays blocky;
- gate tunnelling through a curved stack uses the planar WKB barrier along each normal path (the field is exact, the barrier shape ignores the film's curvature).
## If a solve does not converge
- Sweep in smaller steps, or approach high bias from a converged point.
- Refine the mesh where fields are high (increase Nx/Ny, keep `Graded mesh` on).
- Increase `Gummel its/step` in the Solver page.
- Very high injection (IGBT on-state far above rating) and floating high-voltage regions converge slowly. The status bar says so, and the last converged point is kept.
## Performance
- Use all CPU cores (Solver page).
- 2-D cross-sections: Nz = 3–6 nodes.
- For interactive work keep 3-D meshes under about 200k nodes.
- The triangular mesh has fewer nodes but its multigrid is less effective than the rectangular one, so box geometry runs at a similar speed (faster on large power devices, slower on some others); use it where the geometry needs it.
"""),
("Files and export", """
# Files and export
- File → Save device / Open device: the full editor state (regions, contacts, nets, sheets, mesh, solver and sweep settings) as JSON.
- File → Export solution CSV: every node's x, y, z, φ, n, p, Ec, Ev, J and mobility.
- File → Export sweep CSV: the voltage, the floating-node voltage and the current and resolution floor of every terminal.
- File → Export plot: the visible tab as PNG, PDF or SVG.
- `validate.py`: physics regression battery (python3 validate.py; --templates also solves every template).
- `semimesh.py`: the triangular-prism mesher; it must sit next to main.py.
- `build.sh`: rebuilds the C core on your machine (-march=native).
"""),
("Keyboard shortcuts", """
# Keyboard shortcuts
  F1            Help
  F5            Run (solve at the contact voltages)
  F6            Sweep
  Esc           Stop the running solve or sweep
  Ctrl+E        Show / hide the editor panel
  Ctrl+O        Open device (JSON)
  Ctrl+S        Save device (JSON)
  Ctrl +/-/0    Text size larger / smaller / reset
  Ctrl+Q        Quit
  Mouse wheel   Scroll the panel under the pointer
  Shift+wheel   Scroll sideways
"""),
]

def _clip_poly(P,u0,u1,v0,v1):
    """Polygon P (n,2) clipped to the box [u0,u1] x [v0,v1] (Sutherland-Hodgman)."""
    for axis,val,ge in ((0,u0,True),(0,u1,False),(1,v0,True),(1,v1,False)):
        if len(P)==0: break
        out=[]; n=len(P)
        for i in range(n):
            a,b=P[i],P[(i+1)%n]
            ia=a[axis]>=val if ge else a[axis]<=val; ib=b[axis]>=val if ge else b[axis]<=val
            if ia: out.append(a)
            if ia!=ib: out.append(a+(val-a[axis])/(b[axis]-a[axis])*(b-a))
        P=np.array(out) if out else np.zeros((0,2))
    return P

def _section_polys(sh,box,H,Vv,N,pos):
    """Pieces of (box ∩ outline) - or (box minus outline) for a hole - cut by
    the section plane N = pos, as [(polygon (n,2) in plane coords, holes)]."""
    ax=sh.axis; u,v=[a for a in "xyz" if a!=ax]
    (u0,u1),(v0,v1)=box[u],box[v]
    if N==ax:                                         # section normal to the outline axis
        pl=sh.polyline(2e-3)
        inner=_clip_poly(pl,u0,u1,v0,v1)
        rect=np.array([(u0,v0),(u1,v0),(u1,v1),(u0,v1)])
        tohv=lambda Q:np.c_[Q[:,0] if H==u else Q[:,1], Q[:,1] if Vv==v else Q[:,0]]
        if sh.hole: return [(tohv(rect),[tohv(inner)] if len(inner)>2 else [])]
        return [(tohv(inner),[])] if len(inner)>2 else []
    # the plane contains the outline axis: the outline cut by the line N = pos
    ci=0 if N==u else 1
    pl=sh.polyline(2e-3); c=pl[:,ci]-pos; c2=np.roll(c,-1); q=pl[:,1-ci]; q2=np.roll(q,-1)
    hit=(c*c2<=0)&(c!=c2)
    xs=q[hit]+c[hit]/(c[hit]-c2[hit])*(q2[hit]-q[hit])
    wlo,whi=(box[v] if ci==0 else box[u])
    ivs=[]
    if len(xs)>=2:
        a,b=max(xs.min(),wlo),min(xs.max(),whi)
        ivs=[(wlo,a),(b,whi)] if sh.hole else ([(a,b)] if b>a else [])
    elif sh.hole: ivs=[(wlo,whi)]
    wax=v if ci==0 else u
    out=[]
    for a,b in ivs:
        if b<=a: continue
        (h0,h1)=box[H] if H!=wax else (a,b); (g0,g1)=box[Vv] if Vv!=wax else (a,b)
        out.append((np.array([(h0,g0),(h1,g0),(h1,g1),(h0,g1)]),[]))
    return out

def _poly_patch(pg,**kw):
    """PathPatch of a polygon with holes: pg = (outer (n,2), [inner, ...])."""
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    outer,holes=pg
    verts=[outer,*[h[::-1] for h in holes]]
    codes=[]
    for vv in verts: codes+= [Path.MOVETO]+[Path.LINETO]*(len(vv)-1)
    return PathPatch(Path(np.concatenate(verts),codes),**kw)

def _prism_faces(sh,box,caps=True,nmax=64):
    """3-D faces (plot coordinates x, z, -y) of an outline extruded over its box."""
    ax=sh.axis; u,v=[a for a in "xyz" if a!=ax]
    pl=_clip_poly(sh.polyline(2e-3),box[u][0],box[u][1],box[v][0],box[v][1])
    if len(pl)<3: return []
    if len(pl)>nmax: pl=pl[np.linspace(0,len(pl)-1,nmax).astype(int)]
    a0,a1=box[ax]
    def P(a,q):
        d={ax:a,u:q[0],v:q[1]}; return (d["x"],d["z"],-d["y"])
    faces=[[P(a0,pl[i]),P(a0,pl[(i+1)%len(pl)]),P(a1,pl[(i+1)%len(pl)]),P(a1,pl[i])] for i in range(len(pl))]
    if caps: faces+=[[P(a0,q) for q in pl],[P(a1,q) for q in pl]]
    return faces

def con_extent(c,Lx,Ly,Lz):
    """Physical extent (nm) of a contact: {axis: (lo, hi)}; lo == hi on a face."""
    L={"x":Lx,"y":Ly,"z":Lz}
    if c.face==6:
        return {"x":(c.i0p*Lx,c.i1p*Lx),"y":(c.j0p*Ly,c.j1p*Ly),"z":(c.k0p*Lz,c.k1p*Lz)}
    n="xyz"[c.face//2]; v=0. if c.face%2==0 else L[n]
    a1,a2={"x":("y","z"),"y":("x","z"),"z":("x","y")}[n]
    return {n:(v,v),a1:(c.i0p*L[a1],c.i1p*L[a1]),a2:(c.j0p*L[a2],c.j1p*L[a2])}

def _ram_bytes():
    """Available physical memory (bytes); a large number if unknown."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"): return int(line.split()[1])*1024
    except Exception: pass
    try:
        if sys.platform.startswith("win"):
            class MS(ctypes.Structure):
                _fields_=[("dwLength",ctypes.c_ulong),("dwMemoryLoad",ctypes.c_ulong),
                          ("ullTotalPhys",ctypes.c_ulonglong),("ullAvailPhys",ctypes.c_ulonglong),
                          ("ullTotalPageFile",ctypes.c_ulonglong),("ullAvailPageFile",ctypes.c_ulonglong),
                          ("ullTotalVirtual",ctypes.c_ulonglong),("ullAvailVirtual",ctypes.c_ulonglong),
                          ("ullAvailExtendedVirtual",ctypes.c_ulonglong)]
            m=MS(); m.dwLength=ctypes.sizeof(MS)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)): return int(m.ullAvailPhys)
    except Exception: pass
    return 1<<62

def _map3(x,y,z):
    """Device (x, y=depth, z) -> 3-D plot axes (x, z, -y): top surface up."""
    return x,z,-np.asarray(y) if not np.isscalar(y) else -y

def _depth_axis(ax3):
    ax3.zaxis.set_major_formatter(FuncFormatter(lambda v,p:f"{(-v)+0.:g}"))


# ══════════════════════════════════════════════════════════════
#  APPLICATION
# ══════════════════════════════════════════════════════════════
class App(tk.Tk):
    TEXT_SCALES=[0.8,0.9,1.0,1.15,1.3,1.5,1.75]
    _TABS=[("Structure","td"),("3D Slicer","ts"),("Bands","tb"),("Current","tc"),
           ("I-V","ti"),("Field & Mobility","tf"),("Results","tr")]

    def __init__(self):
        super().__init__()
        self.title(f"SemiSim 3D {VERSION} — drift-diffusion device simulator")
        self.configure(bg=C["bg"])
        self.cfg=self._load_settings()
        sw,sh=self.winfo_screenwidth(),self.winfo_screenheight()
        W=min(1680,max(800,sw-40)); H=min(1020,max(560,sh-90))
        geo=self.cfg.get("geometry")
        try:                                   # a saved window must fit this screen
            import re
            g=re.match(r"(\d+)x(\d+)\+(-?\d+)\+(-?\d+)",geo or "")
            gw,gh,gx,gy=map(int,g.groups())
            if gw>sw or gh>sh or gx<-50 or gy<-50 or gx>sw-150 or gy>sh-150: geo=None
        except Exception: geo=None
        self.geometry(geo if geo else f"{W}x{H}+{max(0,(sw-W)//2)}+{max(0,(sh-H)//4)}")
        self.minsize(760,500)
        # ── state ──
        self.regs:list[Reg]=[]; self.cons:list[Con]=[]; self.sheets:list[Sheet]=[]
        self.sol=None; self.ivd=None; self.busy=False; self.job=None
        self.tname=""; self.tmeta={}
        self.dense_x=[]; self.dense_y=[]; self.dense_z=[]
        self._c_sig=None; self._c_ok=False; self._grid=None; self._mesh_cache=(None,None)
        self._q=queue.Queue(); self._live=[]; self._live_last=0.; self._poll_id=None
        self._helpwin=None; self._browser=None; self._gmap={}
        # domain / mesh
        self.Lx=tk.DoubleVar(value=400.); self.Ly=tk.DoubleVar(value=200.); self.Lz=tk.DoubleVar(value=200.)
        self.Nx=tk.IntVar(value=36); self.Ny=tk.IntVar(value=20); self.Nz=tk.IntVar(value=20)
        self.use_graded=tk.BooleanVar(value=True)
        self.mesh_type=tk.StringVar(value="rect")     # "rect" tensor mesh | "tri" triangular prisms
        self.prism_ax=tk.StringVar(value="auto")      # prism (extrusion) axis of the triangular mesh
        self._prism_cache=(None,None)
        # solver (v7 meanings: CG its / BiCGSTAB its / Gummel its per bias step /
        # Gummel tolerance (V) / relative terminal-current tolerance)
        self.T=tk.DoubleVar(value=300.); self.max_p=tk.IntVar(value=2000)
        self.max_c=tk.IntVar(value=1000); self.max_g=tk.IntVar(value=500)
        self.tol_p=tk.DoubleVar(value=1e-8); self.tol_c=tk.DoubleVar(value=1e-6)
        self.m_bgn=tk.BooleanVar(value=True); self.m_tau=tk.BooleanVar(value=True)
        self.m_vsat=tk.BooleanVar(value=True); self.reuse_sol=tk.BooleanVar(value=True)
        # nanoscale MOS physics (set by each template, see its "models" key)
        self.m_qc=tk.BooleanVar(value=False); self.m_lb=tk.BooleanVar(value=False)
        self.m_tun=tk.BooleanVar(value=False)
        try: ncpu=len(os.sched_getaffinity(0))
        except Exception: ncpu=os.cpu_count() or 1
        self.ncpu=ncpu; self.threads=tk.IntVar(value=int(os.environ.get("OMP_NUM_THREADS",ncpu)))
        # sweep
        self.ivS=tk.DoubleVar(value=-0.5); self.ivE=tk.DoubleVar(value=0.8); self.ivN=tk.IntVar(value=27)
        self.iv_tgt=tk.StringVar(); self.iv_flt=tk.StringVar(value="(none)"); self.iv_out=tk.StringVar(value="auto")
        self.bv_thr=tk.DoubleVar(value=0.); self.iv_refine=tk.BooleanVar(value=True)
        self.iv_live=tk.BooleanVar(value=True)
        # structure view
        self.sec_plane=tk.StringVar(value="XY"); self.sec_pos=tk.DoubleVar(value=0.)
        self.sec_mesh=tk.BooleanVar(value=True); self.sec_3d=tk.BooleanVar(value=True)
        # slicer
        self.sl_axis=tk.StringVar(value="Z"); self.sl_qty=tk.StringVar(value="log10(n)")
        self.sl_pos=tk.IntVar(value=0)
        self._sl_after_id=self._bd_after_id=self._cur_after_id=self._fi_after_id=None
        self._sec_after_id=None
        # band / current / field cuts
        self.bd_axis=tk.StringVar(value="X"); self.bd_p1=tk.IntVar(value=0); self.bd_p2=tk.IntVar(value=0)
        self.show_ef=tk.BooleanVar(value=True); self.show_efnp=tk.BooleanVar(value=True)
        self.show_carriers=tk.BooleanVar(value=True)
        self.show_regions=tk.BooleanVar(value=True)   # region overlay on solution plots
        self.cur_axis=tk.StringVar(value="X"); self.cur_p1=tk.IntVar(value=0); self.cur_p2=tk.IntVar(value=0)
        self.cur_show_arrows2d=tk.BooleanVar(value=True); self.cur_show_arrows3d=tk.BooleanVar(value=True)
        self.fi_axis=tk.StringVar(value="X"); self.fi_p1=tk.IntVar(value=0); self.fi_p2=tk.IntVar(value=0)
        self.fi_show_field=tk.BooleanVar(value=True); self.fi_show_mob=tk.BooleanVar(value=True)
        # region form
        self.rb_mat=tk.StringVar(value="Si"); self.rb_dt=tk.StringVar(value="p")
        self.rb_dp=tk.StringVar(value="1e16"); self.rb_lbl=tk.StringVar(value="")
        self.rb_x0=tk.DoubleVar(); self.rb_y0=tk.DoubleVar(); self.rb_z0=tk.DoubleVar()
        self.rb_x1=tk.DoubleVar(value=200.); self.rb_y1=tk.DoubleVar(value=200.); self.rb_z1=tk.DoubleVar(value=200.)
        # contact form
        self.cb_face=tk.StringVar(value="Y-Min"); self.cb_bc=tk.StringVar(value="Ohmic")
        self.cb_V=tk.DoubleVar(); self.cb_pm=tk.DoubleVar(value=4.05); self.cb_lbl=tk.StringVar()
        self.cb_net=tk.StringVar()
        self.cb_i0=tk.DoubleVar(value=0.); self.cb_i1=tk.DoubleVar(value=100.)
        self.cb_j0=tk.DoubleVar(value=0.); self.cb_j1=tk.DoubleVar(value=100.)
        self.cb_k0=tk.DoubleVar(value=0.); self.cb_k1=tk.DoubleVar(value=100.)
        self.cb_tox=tk.DoubleVar(value=10.); self.cb_nox=tk.DoubleVar(value=0.)
        self.cb_stack=tk.StringVar(value="")
        # text size
        sc=self.cfg.get("text_scale",1.0)
        self.txt_scale=tk.DoubleVar(value=sc if sc in self.TEXT_SCALES else 1.0)
        # defaults for widgets that do not set their own colours (the
        # matplotlib toolbars pick light icons on a dark button background)
        for k,v in (("*Button.background","#0b1422"),("*Button.foreground","#c9d8ec"),
                    ("*Button.activeBackground","#182436"),("*Checkbutton.background","#0b1422"),
                    ("*Checkbutton.foreground","#c9d8ec"),("*Checkbutton.selectColor","#182436"),
                    ("*Label.background","#020609"),("*Label.foreground","#9fb3c8")):
            self.option_add(k,v,"widgetDefault")
        self._fonts(); self._style(); self._build_ui(); self._bind_keys()
        self._apply_text_scale(redraw=False)
        self.protocol("WM_DELETE_WINDOW",self._quit)
        last=self.cfg.get("template")
        self._load(last if last in TEMPLATES else "PN Diode")
        if not self.cfg.get("drawer",True): self.after(50,self._toggle_drawer)

    # ─── settings file ─────────────────────────────────────────
    def _load_settings(self):
        try:
            with open(SETTINGS_FILE) as f: d=json.load(f)
            return d if isinstance(d,dict) else {}
        except Exception: return {}
    def _save_settings(self):
        try:
            w=None
            if self._drawer_on:
                try: w=self.pw.sashpos(0)
                except Exception: w=None
            d=dict(geometry=self.geometry(),text_scale=self.txt_scale.get(),
                   template=self.tname,drawer=self._drawer_on,drawer_w=w or self._drawer_w)
            with open(SETTINGS_FILE,"w") as f: json.dump(d,f)
        except Exception: pass
    def _quit(self):
        if self.busy and self.job is not None: self.job.request_stop()
        self._save_settings(); self.destroy()

    # ─── fonts & style ─────────────────────────────────────────
    def _fonts(self):
        ui=_first_family(self,["Segoe UI","Ubuntu","DejaVu Sans","Noto Sans","Liberation Sans",
                               "Helvetica Neue","Arial","Helvetica"],"TkDefaultFont")
        mono=_first_family(self,["Consolas","Cascadia Mono","DejaVu Sans Mono","Ubuntu Mono",
                                 "Liberation Mono","Noto Sans Mono","Menlo","Courier New"],"TkFixedFont")
        mk=lambda fam,sz,**kw: tkfont.Font(self,family=fam,size=sz,**kw)
        self.F={"ui":mk(ui,9),"uib":mk(ui,9,weight="bold"),"small":mk(ui,8),
                "head":mk(ui,10,weight="bold"),"big":mk(ui,12,weight="bold"),
                "mono":mk(mono,9),"monos":mk(mono,8),"tab":mk(ui,9)}
        self._fbase={k:f.cget("size") for k,f in self.F.items()}
        self._tkbase={}
        for nm in ("TkDefaultFont","TkTextFont","TkMenuFont","TkHeadingFont","TkFixedFont",
                   "TkTooltipFont","TkCaptionFont","TkSmallCaptionFont","TkIconFont"):
            try: self._tkbase[nm]=tkfont.nametofont(nm).cget("size")
            except Exception: pass

    def _style(self):
        s=ttk.Style(self)
        try: s.theme_use("clam")
        except Exception: pass
        bg,pn,bd,tx,sub,ac=C["bg"],C["panel"],C["border"],C["text"],C["sub"],C["accent"]
        s.configure(".",background=pn,foreground=tx,fieldbackground="#07111e",bordercolor=bd,
                    darkcolor=pn,lightcolor=pn,troughcolor=bg,font=self.F["ui"])
        s.configure("TFrame",background=pn)
        s.configure("TLabel",background=pn,foreground=tx,font=self.F["ui"])
        s.configure("TNotebook",background=bg,borderwidth=0,tabmargins=[2,4,2,0])
        s.configure("TNotebook.Tab",background=pn,foreground=sub,padding=[10,4],font=self.F["tab"])
        s.map("TNotebook.Tab",background=[("selected",bg)],foreground=[("selected",ac)])
        s.configure("Drawer.TNotebook",background=pn)
        s.configure("Treeview",background="#07111e",foreground=tx,fieldbackground="#07111e",font=self.F["monos"])
        s.map("Treeview",background=[("selected","#123a5c")],foreground=[("selected","#ffffff")])
        s.configure("Treeview.Heading",background=pn,foreground=ac,font=self.F["small"],relief="flat")
        s.map("Treeview.Heading",background=[("active",bd)])
        s.configure("TProgressbar",background=ac,troughcolor=bg,bordercolor=bd,lightcolor=ac,darkcolor=ac)
        s.configure("TScale",background=pn,troughcolor=bd)
        for st in ("TCheckbutton","TRadiobutton"):
            s.configure(st,background=pn,foreground=tx,font=self.F["ui"],indicatorbackground="#07111e",
                        indicatorforeground=ac,upperbordercolor=sub,lowerbordercolor=sub)
            s.map(st,background=[("active",pn)],foreground=[("disabled",sub)],
                  indicatorbackground=[("pressed",bd),("disabled",pn),("!disabled","#07111e")],
                  indicatorforeground=[("disabled",sub),("!disabled",ac)])
        s.configure("TCombobox",fieldbackground="#07111e",background=pn,foreground=tx,arrowcolor=tx,
                    selectbackground="#07111e",selectforeground=tx)
        s.map("TCombobox",fieldbackground=[("readonly","#07111e")],foreground=[("readonly",tx)],
              selectbackground=[("readonly","#07111e")],selectforeground=[("readonly",tx)])
        s.configure("TSpinbox",fieldbackground="#07111e",foreground=tx,arrowcolor=tx,background=pn)
        s.configure("TScrollbar",background=pn,troughcolor=bg,bordercolor=bd,arrowcolor=sub,
                    lightcolor=pn,darkcolor=pn)
        s.map("TScrollbar",background=[("active",bd)])
        s.configure("TPanedwindow",background=bd)
        s.configure("Sash",sashthickness=6,gripcount=12)
        s.configure("TMenubutton",background=pn,foreground=tx,font=self.F["uib"],padding=[8,3])
        s.map("TMenubutton",background=[("active",bd)])
        s.configure("TSeparator",background=bd)
        self.option_add("*TCombobox*Listbox.background","#07111e")
        self.option_add("*TCombobox*Listbox.foreground",tx)
        self.option_add("*TCombobox*Listbox.selectBackground",ac)
        self.option_add("*TCombobox*Listbox.font",self.F["ui"])
        # the wheel must scroll pages, not silently change combobox values
        for seq in ("<MouseWheel>","<Button-4>","<Button-5>"):
            try: self.unbind_class("TCombobox",seq)
            except Exception: pass

    def _apply_text_scale(self,redraw=True):
        sc=self.txt_scale.get()
        for k,f in self.F.items():
            b=self._fbase[k]; n=int(round(b*sc)); f.configure(size=n if n!=0 else (1 if b>0 else -1))
        for nm,b in self._tkbase.items():
            try:
                n=int(round(b*sc)); tkfont.nametofont(nm).configure(size=n if n!=0 else (1 if b>0 else -1))
            except Exception: pass
        ls=self.F["monos"].metrics("linespace")
        ttk.Style(self).configure("Treeview",rowheight=int(ls*1.3)+2)
        self.after_idle(self._fit_trees); self.after(120,self._fit_drawer)
        for fig in getattr(self,"_figs",[]):
            self._scale_fig(fig,sc)
        if redraw: self._redraw_all()

    @staticmethod
    def _scale_fig(fig,sc):
        """Plot text follows the UI text size: change the figure dpi while
        keeping its pixel size (the Tk backend re-derives dpi from
        _original_dpi when the screen's pixel ratio changes)."""
        try:
            W,H=fig.bbox.width,fig.bbox.height
            r=getattr(fig.canvas,"_device_pixel_ratio",1.) or 1.
            fig._original_dpi=100.*sc; dpi=100.*sc*r
            fig._set_dpi(dpi,forward=False)
            if W>2 and H>2: fig.set_size_inches(W/dpi,H/dpi,forward=False)
        except Exception: pass

    def _text_step(self,d):
        L=self.TEXT_SCALES; cur=self.txt_scale.get()
        i=min(range(len(L)),key=lambda k:abs(L[k]-cur))
        i=0 if d is None else max(0,min(len(L)-1,i+d))
        self.txt_scale.set(L[i] if d is not None else 1.0); self._apply_text_scale()

    # ─── small widget factories ────────────────────────────────
    def _sec(self,p,t):
        f=tk.LabelFrame(p,text=f" {t} ",bg=C["panel"],fg=C["accent"],font=self.F["uib"],
                        relief="groove",bd=1,labelanchor="nw")
        f.pack(fill=tk.X,padx=6,pady=4,anchor="n"); return f
    def _L(self,p,t,r,c,fg=None,**kw):
        l=tk.Label(p,text=t,bg=C["panel"],fg=fg or C["sub"],font=self.F["ui"])
        l.grid(row=r,column=c,sticky="w",padx=(4,2),pady=2,**kw); return l
    def _E(self,p,v,r,c,w=8,**kw):
        e=tk.Entry(p,textvariable=v,bg="#07111e",fg=C["text"],insertbackground=C["text"],relief="flat",
                   font=self.F["mono"],width=w,highlightbackground=C["border"],highlightcolor=C["accent"],
                   highlightthickness=1)
        e.grid(row=r,column=c,padx=2,pady=2,sticky="w",**kw); return e
    def _CB(self,p,v,opts,r,c,w=10,**kw):
        cb=ttk.Combobox(p,textvariable=v,values=list(opts),state="readonly",width=w,font=self.F["ui"])
        cb.grid(row=r,column=c,padx=2,pady=2,sticky="w",**kw); return cb
    def _btn(self,p,t,cmd,bg=None,fg="white",**kw):
        fn=kw.pop("font",self.F["ui"])
        return tk.Button(p,text=t,command=cmd,bg=bg or C["accent"],fg=fg,activebackground=C["border"],
                         activeforeground="white",relief="flat",cursor="hand2",font=fn,padx=8,pady=3,
                         disabledforeground="#6b7c90",**kw)
    def _hint(self,p,t,**kw):
        """Wrapped help text inside an editor page.  The wrap width follows
        the VISIBLE drawer width (wrapping to the label's own width would let
        a long line widen the page instead)."""
        l=tk.Label(p,text=t,bg=C["panel"],fg=C["sub"],font=self.F["small"],justify="left",anchor="w",
                   wraplength=self._hint_wrap())
        l.pack(fill=tk.X,padx=8,pady=(0,4),**kw)
        self._hints=getattr(self,"_hints",[])+[l]
        return l
    def _hint_wrap(self):
        try: w=self.drawer.winfo_width()
        except Exception: w=0
        return max(220,(w if w>50 else 420)-70)
    def _rewrap(self,_=None):
        wl=self._hint_wrap()
        for l in getattr(self,"_hints",[]):
            try: l.configure(wraplength=wl)
            except Exception: pass
    def _tree(self,parent,cols,widths,height=7,anchor="center",fixed=True):
        """Treeview with both scrollbars in its own frame.  fixed=True: the
        frame does not grow with the table's natural width (a 10-column table
        must not widen the whole editor page) - the table scrolls instead."""
        fr=tk.Frame(parent,bg=C["panel"])
        tv=ttk.Treeview(fr,columns=cols,show="headings",height=height,selectmode="browse")
        for c2,w in zip(cols,widths):
            tv.heading(c2,text=c2); tv.column(c2,width=w,minwidth=30,anchor=anchor,stretch=False)
        vs=ttk.Scrollbar(fr,orient="vertical",command=tv.yview)
        hs=ttk.Scrollbar(fr,orient="horizontal",command=tv.xview)
        tv.configure(yscrollcommand=vs.set,xscrollcommand=hs.set)
        tv.grid(row=0,column=0,sticky="nsew"); vs.grid(row=0,column=1,sticky="ns"); hs.grid(row=1,column=0,sticky="ew")
        fr.rowconfigure(0,weight=1); fr.columnconfigure(0,weight=1)
        if fixed:
            fr.grid_propagate(False)
            self._trees=getattr(self,"_trees",[])+[(fr,tv,hs)]
            self.after_idle(self._fit_trees)
        return fr,tv
    def _fit_trees(self):
        for fr,tv,hs in getattr(self,"_trees",[]):
            try:
                fr.update_idletasks()
                fr.configure(width=240,height=tv.winfo_reqheight()+hs.winfo_reqheight()+2)
            except Exception: pass
    def _mfig(self,p,fs=(11,5.5)):
        fig,cv=mfig(p,fs)
        cv.get_tk_widget().wheel_native=True
        self._scale_fig(fig,self.txt_scale.get())
        self._figs=getattr(self,"_figs",[])+[fig]
        return fig,cv

    # ─── BUILD UI ──────────────────────────────────────────────
    def _build_ui(self):
        self._menus()
        self._toolbar()
        sb=tk.Frame(self,bg="#020508"); sb.pack(fill=tk.X,side=tk.BOTTOM)
        self.sv=tk.StringVar(value="Ready")
        self.sv_lbl=tk.Label(sb,textvariable=self.sv,bg="#020508",fg=C["sub"],font=self.F["small"],
                             anchor="w",padx=8)
        self.sv_lbl.pack(side=tk.LEFT,fill=tk.X,expand=True)
        self.cv2=tk.StringVar(value="")
        tk.Label(sb,textvariable=self.cv2,bg="#020508",fg=C["green"],font=self.F["small"],padx=8).pack(side=tk.RIGHT)
        self.pw=ttk.PanedWindow(self,orient=tk.HORIZONTAL)
        self.pw.pack(fill=tk.BOTH,expand=True,padx=2,pady=(0,2))
        self.drawer=tk.Frame(self.pw,bg=C["panel"],width=440)
        self.plots=tk.Frame(self.pw,bg=C["bg"])
        self.pw.add(self.drawer,weight=0); self.pw.add(self.plots,weight=1)
        self._drawer_on=True; self._drawer_w=None
        self._build_drawer(self.drawer); self._right(self.plots)
        self.drawer.bind("<Configure>",lambda e:self._debounce("_wrap_after_id",self._rewrap,60),add="+")
        self._sash_init=False
        def first_layout(e):
            if not self._sash_init and e.width>300:
                self._sash_init=True
                w=int(self.cfg.get("drawer_w") or min(480,max(340,0.32*e.width)))
                self.after_idle(lambda:self._set_sash(w))
        self.pw.bind("<Configure>",first_layout,add="+")

    def _fit_drawer(self):
        """After a text-size change: widen the drawer so its tabs fit
        (at most 55 % of the window)."""
        if not self._drawer_on: return
        try:
            need=self.dnb.winfo_reqwidth()+6; cur=self.pw.sashpos(0)
            mx=int(0.55*self.pw.winfo_width())
            if cur<need: self._set_sash(min(need,mx))
        except Exception: pass
    def _set_sash(self,x):
        try: self.pw.sashpos(0,x)
        except Exception: pass

    def _menus(self):
        mkw=dict(tearoff=0,bg=C["panel"],fg=C["text"],activebackground=C["accent"],
                 activeforeground="white",font=self.F["ui"])
        mb=tk.Menu(self,**mkw); self.config(menu=mb)
        fm=tk.Menu(mb,**mkw)
        fm.add_command(label="Open device…",accelerator="Ctrl+O",command=self._open_device)
        fm.add_command(label="Save device as…",accelerator="Ctrl+S",command=self._save_device)
        fm.add_separator()
        fm.add_command(label="Export solution CSV…",command=self._xcsv)
        fm.add_command(label="Export sweep CSV…",command=self._xiv_csv)
        fm.add_command(label="Export visible plot…",command=self._xplt)
        fm.add_separator(); fm.add_command(label="Quit",accelerator="Ctrl+Q",command=self._quit)
        mb.add_cascade(label="File",menu=fm)
        dm=tk.Menu(mb,**mkw); self._fill_template_menu(dm,mkw)
        mb.add_cascade(label="Devices",menu=dm)
        rm=tk.Menu(mb,**mkw)
        rm.add_command(label="▶ Run",accelerator="F5",command=self._run)
        rm.add_command(label="↔ Sweep",accelerator="F6",command=self._sweep)
        rm.add_command(label="■ Stop",accelerator="Esc",command=self._stop)
        rm.add_separator()
        rm.add_checkbutton(label="Continue from previous solution",variable=self.reuse_sol)
        rm.add_command(label="Forget solution (next run starts from equilibrium)",command=self._forget)
        mb.add_cascade(label="Run",menu=rm)
        vm=tk.Menu(mb,**mkw)
        vm.add_command(label="Show / hide editor panel",accelerator="Ctrl+E",command=self._toggle_drawer)
        ts=tk.Menu(vm,**mkw)
        for sc in self.TEXT_SCALES:
            ts.add_radiobutton(label=f"{int(sc*100)} %",variable=self.txt_scale,value=sc,
                               command=self._apply_text_scale)
        vm.add_cascade(label="Text size",menu=ts)
        vm.add_command(label="Larger text",accelerator="Ctrl++",command=lambda:self._text_step(1))
        vm.add_command(label="Smaller text",accelerator="Ctrl+-",command=lambda:self._text_step(-1))
        vm.add_separator()
        for i,(nm,_) in enumerate(self._TABS):
            vm.add_command(label=f"{nm.strip()} tab",command=lambda i=i:self.nb.select(i))
        mb.add_cascade(label="View",menu=vm)
        hm=tk.Menu(mb,**mkw)
        hm.add_command(label="Help contents",accelerator="F1",command=self._help)
        hm.add_command(label="Keyboard shortcuts",command=lambda:self._help("Keyboard shortcuts"))
        hm.add_command(label="Template library…",command=self._browse_templates)
        hm.add_separator(); hm.add_command(label="About SemiSim 3D",command=self._about)
        mb.add_cascade(label="Help",menu=hm)

    def _fill_template_menu(self,m,mkw):
        for cat,names in template_groups():
            sub=tk.Menu(m,**mkw)
            for n in names: sub.add_command(label=n,command=lambda n=n:self._load(n))
            m.add_cascade(label=cat,menu=sub)
        m.add_separator(); m.add_command(label="Browse templates with descriptions…",command=self._browse_templates)

    def _toolbar(self):
        tb=tk.Frame(self,bg=C["panel"],highlightthickness=1,highlightbackground=C["border"])
        tb.pack(fill=tk.X,side=tk.TOP)
        self.toolbar=tb
        b=self._btn(tb,"☰",self._toggle_drawer,bg=C["border"],font=self.F["big"]); b.pack(side=tk.LEFT,padx=(4,2),pady=3)
        Tip(b,"Show / hide the editor panel (Ctrl+E)",self.F["small"])
        tk.Label(tb,text="Device",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(8,2))
        self.tmpl_btn=ttk.Menubutton(tb,text="—",style="TMenubutton",width=26)
        mkw=dict(tearoff=0,bg=C["panel"],fg=C["text"],activebackground=C["accent"],
                 activeforeground="white",font=self.F["ui"])
        tm=tk.Menu(self.tmpl_btn,**mkw); self._fill_template_menu(tm,mkw)
        self.tmpl_btn["menu"]=tm; self.tmpl_btn.pack(side=tk.LEFT,padx=2,pady=3)
        Tip(self.tmpl_btn,"Template library by category",self.F["small"])
        ttk.Separator(tb,orient="vertical").pack(side=tk.LEFT,fill=tk.Y,padx=6,pady=4)
        self.b_run=self._btn(tb,"▶ Run",self._run,bg="#0a8f5a",font=self.F["uib"])
        self.b_run.pack(side=tk.LEFT,padx=2,pady=3)
        Tip(self.b_run,"Solve at the contact voltages (F5)",self.F["small"])
        self.b_swp=self._btn(tb,"↔ Sweep",self._sweep,bg="#7a3fc4",font=self.F["uib"])
        self.b_swp.pack(side=tk.LEFT,padx=2,pady=3)
        Tip(self.b_swp,"Run the sweep defined in the Sweep page (F6)",self.F["small"])
        self.b_stop=self._btn(tb,"■ Stop",self._stop,bg="#a3243b",font=self.F["uib"],state="disabled")
        self.b_stop.pack(side=tk.LEFT,padx=2,pady=3)
        Tip(self.b_stop,"Stop the solver (Esc) - the last converged bias point is kept",self.F["small"])
        hb=self._btn(tb,"?",self._help,bg=C["border"],font=self.F["uib"]); hb.pack(side=tk.RIGHT,padx=4,pady=3)
        Tip(hb,"Help (F1)",self.F["small"])
        self.pg=ttk.Progressbar(tb,mode="determinate",maximum=100.,length=150)
        self.pg.pack(side=tk.LEFT,padx=(10,4),pady=3)
        self.pg_txt=tk.StringVar(value="")
        tk.Label(tb,textvariable=self.pg_txt,bg=C["panel"],fg=C["yellow"],font=self.F["small"],
                 anchor="w").pack(side=tk.LEFT,fill=tk.X,expand=True,padx=4)

    def _bind_keys(self):
        self.bind_all("<F1>",lambda e:self._help())
        self.bind_all("<F5>",lambda e:self._run())
        self.bind_all("<F6>",lambda e:self._sweep())
        self.bind_all("<Escape>",lambda e:self._stop())
        for seq,fn in [("<Control-e>",self._toggle_drawer),("<Control-E>",self._toggle_drawer),
                       ("<Control-o>",self._open_device),("<Control-s>",self._save_device),
                       ("<Control-q>",self._quit)]:
            self.bind_all(seq,lambda e,fn=fn:(fn(),"break")[1])
        for seq in ("<Control-plus>","<Control-equal>","<Control-KP_Add>"):
            self.bind_all(seq,lambda e:self._text_step(1))
        for seq in ("<Control-minus>","<Control-KP_Subtract>"):
            self.bind_all(seq,lambda e:self._text_step(-1))
        self.bind_all("<Control-0>",lambda e:self._text_step(None))
        for seq in ("<MouseWheel>","<Shift-MouseWheel>","<Button-4>","<Button-5>",
                    "<Shift-Button-4>","<Shift-Button-5>"):
            self.bind_all(seq,self._wheel,add="+")

    def _wheel(self,e):
        """Scroll the ScrollFrame under the pointer (Windows/macOS MouseWheel
        and X11 buttons 4/5; Shift = sideways).  Widgets that scroll by
        themselves (tables, text, plots, sliders) are left alone."""
        try: w=self.winfo_containing(e.x_root,e.y_root)
        except Exception: return
        if w is None: return
        if getattr(e,"num",None)==4: step=-1
        elif getattr(e,"num",None)==5: step=1
        else:
            d=getattr(e,"delta",0)
            if not d: return
            step=-max(1,abs(d)//120) if d>0 else max(1,abs(d)//120)
        horiz=bool(e.state & 0x0001)
        x=w
        while x is not None:
            if getattr(x,"wheel_native",False): return
            if isinstance(x,(tk.Text,tk.Listbox,ttk.Treeview,ttk.Scale,tk.Scale)): return
            if getattr(x,"wheel_target",False):
                (x.xview_scroll if horiz else x.yview_scroll)(step,"units"); return "break"
            x=getattr(x,"master",None)

    def _toggle_drawer(self):
        if self._drawer_on:
            try: self._drawer_w=self.pw.sashpos(0)
            except Exception: self._drawer_w=None
            self.pw.forget(self.drawer); self._drawer_on=False
        else:
            self.pw.insert(0,self.drawer,weight=0); self._drawer_on=True
            w=self._drawer_w or 440
            self.after(30,lambda:self._set_sash(w))

    # ─── DRAWER (editor pages) ─────────────────────────────────
    def _build_drawer(self,p):
        nb=ttk.Notebook(p,style="Drawer.TNotebook"); nb.pack(fill=tk.BOTH,expand=True)
        self.dnb=nb
        pages={}
        for key,title in [("dev","Device"),("reg","Regions"),("con","Contacts"),
                          ("sol","Solver"),("swp","Sweep")]:
            sf=ScrollFrame(nb,C["panel"]); nb.add(sf,text=title); pages[key]=sf.inner
        self._page_device(pages["dev"]); self._page_regions(pages["reg"])
        self._page_contacts(pages["con"]); self._page_solver(pages["sol"]); self._page_sweep(pages["swp"])

    def _page_device(self,p):
        fr=self._sec(p,"Template")
        self.t_name=tk.Label(fr,text="—",bg=C["panel"],fg=C["accent"],font=self.F["head"],anchor="w",justify="left")
        self.t_name.pack(fill=tk.X,padx=6,pady=(2,0))
        self.t_cat=tk.Label(fr,text="",bg=C["panel"],fg=C["sub"],font=self.F["small"],anchor="w")
        self.t_cat.pack(fill=tk.X,padx=6)
        tf=tk.Frame(fr,bg=C["panel"]); tf.pack(fill=tk.BOTH,expand=True,padx=4,pady=4)
        self.t_doc=tk.Text(tf,height=9,width=30,wrap="word",bg="#07111e",fg=C["text"],relief="flat",
                           font=self.F["small"],padx=6,pady=4,highlightthickness=0)
        ds=ttk.Scrollbar(tf,orient="vertical",command=self.t_doc.yview); self.t_doc.configure(yscrollcommand=ds.set)
        self.t_doc.pack(side=tk.LEFT,fill=tk.BOTH,expand=True); ds.pack(side=tk.RIGHT,fill=tk.Y)
        bf=tk.Frame(fr,bg=C["panel"]); bf.pack(fill=tk.X,padx=4,pady=(0,4))
        self._btn(bf,"Browse templates…",self._browse_templates,bg=C["border"]).pack(side=tk.LEFT,padx=2)
        self._btn(bf,"Reload template",lambda:self._load(self.tname),bg=C["border"]).pack(side=tk.LEFT,padx=2)

        fr=self._sec(p,"Domain & mesh")
        fm=tk.Frame(fr,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"Size (nm)",0,0)
        for k,(t,v) in enumerate([("Lx",self.Lx),("Ly",self.Ly),("Lz",self.Lz)]):
            self._L(fm,t,0,1+2*k); self._E(fm,v,0,2+2*k,8)
        self._L(fm,"Base nodes",1,0)
        for k,(t,v) in enumerate([("Nx",self.Nx),("Ny",self.Ny),("Nz",self.Nz)]):
            self._L(fm,t,1,1+2*k); self._E(fm,v,1,2+2*k,8)
        ttk.Checkbutton(fm,text="Graded mesh (dense at junctions, interfaces, gates)",
                        variable=self.use_graded).grid(row=2,column=0,columnspan=7,sticky="w",padx=4,pady=3)
        self._L(fm,"Mesh",3,0)
        for k,(txt,val) in enumerate((("Rectangular","rect"),("Triangular (prisms)","tri"))):
            ttk.Radiobutton(fm,text=txt,variable=self.mesh_type,value=val,
                            command=self._mesh_changed).grid(row=3,column=1+2*k,columnspan=2,sticky="w",padx=2)
        self._L(fm,"Prism axis",4,0)
        cb=self._CB(fm,self.prism_ax,["auto","x","y","z"],4,1,w=6,columnspan=2)
        cb.bind("<<ComboboxSelected>>",lambda e:self._mesh_changed())
        bf=tk.Frame(fr,bg=C["panel"]); bf.pack(fill=tk.X,padx=4)
        self._btn(bf,"Update mesh & structure",self._mesh_update,bg=C["accent"]).pack(side=tk.LEFT,padx=2,pady=2)
        self.mesh_lbl=tk.Label(fr,text="",bg=C["panel"],fg=C["green"],font=self.F["mono"],anchor="w",justify="left")
        self.mesh_lbl.pack(fill=tk.X,padx=8,pady=2)
        self.dense_lbl=self._hint(fr,"")
        self._hint(fr,"y is depth (y = 0 is the top surface). For a 2-D cross-section keep Lz small in "
                      "nodes (Nz = 3-6); currents are then per Lz of width.")
        self._hint(fr,"Triangular (prisms): the cross-section normal to the prism axis is a Delaunay "
                      "triangulation, extruded along the axis. Round and sloped outlines (the shaped "
                      "templates: GAA nanowire, rounded FinFET) are meshed as they are, not as staircases, "
                      "and fine rows stay near the interfaces that need them. Box geometry gives the same "
                      "results as the rectangular mesh with fewer nodes. The graded setting does not apply.")
        fr=self._sec(p,"Interface sheet charges")
        self.sheet_lbl=self._hint(fr,"none")

    def _page_regions(self,p):
        fr=self._sec(p,"Regions (applied in order - later rows override)")
        cols=("#","Label","Mat","Type","Doping cm⁻³","x0","x1","y0","y1","z0","z1")
        tf,self.rtree=self._tree(fr,cols,[30,170,56,40,78,52,52,52,52,52,52],height=9)
        tf.pack(fill=tk.BOTH,expand=True,padx=4,pady=4)
        self.rtree.bind("<<TreeviewSelect>>",lambda _:self._loadreg())
        fr2=self._sec(p,"Edit region")
        fm=tk.Frame(fr2,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"Material",0,0); self._CB(fm,self.rb_mat,list(MATS),0,1,w=8)
        self._L(fm,"Type",0,2); self._CB(fm,self.rb_dt,["n","p","i"],0,3,w=3)
        self._L(fm,"Doping (cm⁻³)",1,0); self._E(fm,self.rb_dp,1,1,10)
        self._L(fm,"x0 / x1",2,0); self._E(fm,self.rb_x0,2,1,8); self._E(fm,self.rb_x1,2,2,8,columnspan=2)
        self._L(fm,"y0 / y1",3,0); self._E(fm,self.rb_y0,3,1,8); self._E(fm,self.rb_y1,3,2,8,columnspan=2)
        self._L(fm,"z0 / z1",4,0); self._E(fm,self.rb_z0,4,1,8); self._E(fm,self.rb_z1,4,2,8,columnspan=2)
        self._L(fm,"Label",5,0); self._E(fm,self.rb_lbl,5,1,26,columnspan=3)
        bf=tk.Frame(fr2,bg=C["panel"]); bf.pack(fill=tk.X,padx=4,pady=4)
        items=[("+ Add",self._addreg,"#0a8f5a"),("✎ Update",self._editreg,C["accent"]),
               ("⧉ Duplicate",self._dupreg,"#0e7490"),("✕ Delete",self._delreg,"#a3243b"),
               ("▲ Up",self._regup,C["border"]),("▼ Down",self._regdn,C["border"]),
               ("Fill domain",self._fillreg,"#4f46e5"),("Clear all",self._clrreg,"#444")]
        for k,(t2,cmd,bg) in enumerate(items):
            self._btn(bf,t2,cmd,bg=bg).grid(row=k//4,column=k%4,padx=2,pady=2,sticky="ew")
        for c in range(4): bf.columnconfigure(c,weight=1)
        self._hint(fr2,"Select a row to load it into the form, edit, then Update. Coordinates in nm; "
                       "'i' = intrinsic. Insulators (SiO2, Al2O3) ignore the doping.")

    def _page_contacts(self,p):
        fr=self._sec(p,"Contacts")
        cols=("#","Label","Net","Face","BC","V","φm","EOT (stack)","Nox","Range (nm)")
        tf,self.ctree=self._tree(fr,cols,[30,150,60,56,60,56,46,90,52,260],height=8)
        tf.pack(fill=tk.BOTH,expand=True,padx=4,pady=4)
        self.ctree.bind("<<TreeviewSelect>>",lambda _:self._loadcon())
        fr2=self._sec(p,"Edit contact")
        fm=tk.Frame(fr2,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"Face",0,0); self._CB(fm,self.cb_face,Con.FNAME,0,1,w=8)
        self._L(fm,"BC",0,2); self._CB(fm,self.cb_bc,Con.BCNAME,0,3,w=9)
        self._L(fm,"V (V)",1,0); self._E(fm,self.cb_V,1,1,9)
        self._L(fm,"φm (eV)",1,2); self._E(fm,self.cb_pm,1,3,7)
        self._L(fm,"Net",2,0); self._E(fm,self.cb_net,2,1,9)
        self._L(fm,"(same net = wired together)",2,2,columnspan=2)
        self.con_axis_lbl=tk.Label(fm,text="",bg=C["panel"],fg=C["yellow"],font=self.F["uib"])
        self.con_axis_lbl.grid(row=3,column=0,columnspan=4,sticky="w",padx=4,pady=(4,0))
        self.cb_face.trace_add("write",self._upd_con_axes)
        self.ci_lbl=self._L(fm,"i0 / i1 %",4,0); self._E(fm,self.cb_i0,4,1,8); self._E(fm,self.cb_i1,4,2,8)
        self.cj_lbl=self._L(fm,"j0 / j1 %",5,0); self._E(fm,self.cb_j0,5,1,8); self._E(fm,self.cb_j1,5,2,8)
        self.ck_lbl=self._L(fm,"k0 / k1 %",6,0); self._E(fm,self.cb_k0,6,1,8); self._E(fm,self.cb_k1,6,2,8)
        self._L(fm,"Gate tox (nm)",7,0); self._E(fm,self.cb_tox,7,1,8)
        self._L(fm,"Nox (q/cm²)",7,2); self._E(fm,self.cb_nox,7,3,9)
        self._L(fm,"Gate stack",8,0); self._E(fm,self.cb_stack,8,1,26,columnspan=3)
        self._L(fm,"Label",9,0); self._E(fm,self.cb_lbl,9,1,26,columnspan=3)
        bf=tk.Frame(fr2,bg=C["panel"]); bf.pack(fill=tk.X,padx=4,pady=4)
        items=[("+ Add",self._addcon,"#0a8f5a"),("✎ Update",self._editcon,C["accent"]),
               ("⧉ Duplicate",self._dupcon,"#0e7490"),("✕ Delete",self._delcon,"#a3243b"),
               ("Full face",self._full_face,C["border"]),("Clear all",self._clrcon,"#444")]
        for k,(t2,cmd,bg) in enumerate(items):
            self._btn(bf,t2,cmd,bg=bg).grid(row=k//3,column=k%3,padx=2,pady=2,sticky="ew")
        for c in range(3): bf.columnconfigure(c,weight=1)
        self._hint(fr2,"Ranges are % of the face's two in-plane axes (Box: % of X, Y and Z). "
                       "Ohmic = charge-neutral contact; Schottky uses φm; Gate on a semiconductor face = "
                       "oxide of thickness tox (EOT) with fixed charge Nox; Gate on oxide nodes or as a Box = "
                       "metal electrode. Gate stack (face gates, optional): layers from the semiconductor "
                       "outwards in nm, e.g. 'SiO2 0.5, HfO2 2.5' - it replaces tox (EOT computed from the "
                       "permittivities) and is what gate tunnelling sees; empty = SiO2 of thickness tox. "
                       "Contacts sharing a Net always have the same voltage - editing one member's V "
                       "updates the whole net.")
        self._upd_con_axes()

    def _page_solver(self,p):
        fr=self._sec(p,"Physics")
        fm=tk.Frame(fr,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"Temperature (K)",0,0); self._E(fm,self.T,0,1,8)
        for k,(txt,var) in enumerate([("Band-gap narrowing (Klaassen)",self.m_bgn),
                                      ("SRH lifetime τ(N)",self.m_tau),
                                      ("Velocity saturation",self.m_vsat),
                                      ("Quantum confinement at oxide interfaces (MLDA)",self.m_qc),
                                      ("Surface mobility (Lombardi, Si)",self.m_lb),
                                      ("Gate tunnelling (direct + Fowler-Nordheim)",self.m_tun)]):
            ttk.Checkbutton(fm,text=txt,variable=var).grid(row=1+k,column=0,columnspan=3,sticky="w",padx=4,pady=1)
        self._hint(fr,"The last three matter for thin-oxide MOS devices: MLDA pushes the inversion layer "
                      "~1 nm off the oxide (V_th up, gate capacitance down), Lombardi lowers the channel "
                      "mobility with the normal field, and tunnelling lets current flow through the gate "
                      "stack (see each gate's stack on the Contacts page). Templates set these switches.")
        fr=self._sec(p,"Iterations & tolerances")
        fm=tk.Frame(fr,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        rows=[("Gummel iterations per bias step",self.max_g),("Gummel tolerance (V)",self.tol_p),
              ("Terminal-current rel. tolerance",self.tol_c),("CG iterations (Poisson)",self.max_p),
              ("BiCGSTAB iterations (continuity)",self.max_c)]
        for k,(t,v) in enumerate(rows):
            self._L(fm,t,k,0); self._E(fm,v,k,1,10)
        fr=self._sec(p,"Performance")
        fm=tk.Frame(fr,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"CPU threads",0,0)
        sp=ttk.Spinbox(fm,from_=1,to=max(1,self.ncpu),textvariable=self.threads,width=5,
                       command=self._set_threads,font=self.F["mono"])
        sp.grid(row=0,column=1,padx=2,pady=2,sticky="w"); sp.bind("<Return>",lambda e:self._set_threads())
        self._L(fm,f"of {self.ncpu}",0,2)
        ttk.Checkbutton(fm,text="Continue from the previous solution when only voltages changed",
                        variable=self.reuse_sol).grid(row=1,column=0,columnspan=3,sticky="w",padx=4,pady=2)
        self._hint(fr,"Newton-Poisson (multigrid CG) + Scharfetter-Gummel continuity (multigrid BiCGSTAB), "
                      "Anderson-accelerated Gummel map, adaptive bias continuation. Defaults suit every "
                      "template; raise 'Gummel iterations' if a hard bias point stops early.")

    def _page_sweep(self,p):
        fr=self._sec(p,"Sweep")
        fm=tk.Frame(fr,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"Swept terminal",0,0)
        self.cb_tgt=self._CB(fm,self.iv_tgt,[],0,1,w=26,columnspan=3)
        self._L(fm,"V start",1,0); self._E(fm,self.ivS,1,1,9)
        self._L(fm,"V end",1,2); self._E(fm,self.ivE,1,3,9)
        self._L(fm,"Points",2,0); self._E(fm,self.ivN,2,1,9)
        self._L(fm,"Floating (I = 0)",3,0)
        self.cb_flt=self._CB(fm,self.iv_flt,["(none)"],3,1,w=26,columnspan=3)
        self._L(fm,"Plot current of",4,0)
        self.cb_out=self._CB(fm,self.iv_out,["auto"],4,1,w=26,columnspan=3)
        self._L(fm,"BV knee (A/cm²)",5,0); self._E(fm,self.bv_thr,5,1,9)
        self._L(fm,"0 = off",5,2)
        ttk.Checkbutton(fm,text="Adaptive refinement of steep transfer curves",
                        variable=self.iv_refine).grid(row=6,column=0,columnspan=4,sticky="w",padx=4,pady=1)
        ttk.Checkbutton(fm,text="Live plot while sweeping",
                        variable=self.iv_live).grid(row=7,column=0,columnspan=4,sticky="w",padx=4,pady=1)
        bf=tk.Frame(fr,bg=C["panel"]); bf.pack(fill=tk.X,padx=4,pady=4)
        self._btn(bf,"↔ Run sweep (F6)",self._sweep,bg="#7a3fc4",font=self.F["uib"]).pack(side=tk.LEFT,padx=2)
        self._hint(fr,"The swept terminal may be a whole net (all members move together). A floating "
                      "terminal is solved for zero current at every point - e.g. an inverter output - which "
                      "turns the sweep into a transfer curve. All other contacts stay at the voltages in the "
                      "Contacts page.")


    # ─── PLOT TABS ─────────────────────────────────────────────
    def _right(self,p):
        nb=ttk.Notebook(p); nb.pack(fill=tk.BOTH,expand=True); self.nb=nb
        for name,attr in self._TABS:
            f=tk.Frame(nb,bg=C["bg"]); nb.add(f,text=f" {name} "); setattr(self,attr,f)
        self._init_device_tab(); self._init_slicer_tab(); self._init_band_tab()
        self._init_curr_tab(); self._init_fi_tab(); self._init_iv_tab(); self._init_results_tab()

    def _strip(self,parent):
        f=tk.Frame(parent,bg=C["panel"]); f.pack(fill=tk.X,padx=4,pady=(3,1)); return f
    def _radio_axes(self,parent,var,cmd,label="Cut axis",vals=("X","Y","Z")):
        tk.Label(parent,text=label,bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(4,2))
        for v in vals:
            ttk.Radiobutton(parent,text=v,variable=var,value=v,command=cmd).pack(side=tk.LEFT,padx=2)
    def _check(self,parent,text,var,cmd,fg=None):
        tk.Checkbutton(parent,text=text,variable=var,command=cmd,bg=C["panel"],fg=fg or C["text"],
                       selectcolor="#07111e",activebackground=C["panel"],activeforeground=fg or C["text"],
                       font=self.F["ui"],highlightthickness=0).pack(side=tk.LEFT,padx=4)

    def _init_device_tab(self):
        cf=self._strip(self.td)
        self._radio_axes(cf,self.sec_plane,self._sec_changed,"Cross-section",("XY","YZ","XZ"))
        self.sec_lbl=tk.Label(cf,text="",bg=C["panel"],fg=C["yellow"],font=self.F["mono"],width=13,anchor="w")
        self.sec_lbl.pack(side=tk.LEFT,padx=4)
        self.sec_scale=ttk.Scale(cf,from_=0,to=100,orient=tk.HORIZONTAL,variable=self.sec_pos,
                                 command=lambda _:self._sec_moved())
        self.sec_scale.pack(side=tk.LEFT,fill=tk.X,expand=True,padx=4)
        self._check(cf,"mesh",self.sec_mesh,self._ddev)
        self._check(cf,"3-D view",self.sec_3d,self._ddev)
        self.fd,self.cvd=self._mfig(self.td,fs=(13,6))

    def _sec_changed(self):
        self._sec_range(None); self._ddev()
    def _sec_range(self,pos=None):
        L={"XY":self.Lz,"YZ":self.Lx,"XZ":self.Ly}[self.sec_plane.get()]
        try: Lm=float(L.get())
        except Exception: Lm=1.
        self.sec_scale.configure(from_=0.,to=Lm)
        self.sec_pos.set(Lm/2. if pos is None else min(max(pos,0.),Lm))
        self._sec_lbl_upd()
    def _sec_lbl_upd(self):
        n={"XY":"z","YZ":"x","XZ":"y"}[self.sec_plane.get()]
        self.sec_lbl.config(text=f"{n} = {self.sec_pos.get():.1f} nm")
    def _sec_moved(self):
        self._sec_lbl_upd(); self._debounce("_sec_after_id",self._ddev,180)

    def _init_slicer_tab(self):
        cf=self._strip(self.ts)
        self._radio_axes(cf,self.sl_axis,self._upd_slice,"Slice normal")
        tk.Label(cf,text="Quantity",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(12,2))
        qtys=["phi","log10(n)","log10(p)","Ec","Ev","Efn","Efp","log10|J|","mobility µn","log10|R|","net doping",
              "log10|J gate|","Λn quantum"]
        cb=ttk.Combobox(cf,textvariable=self.sl_qty,values=qtys,state="readonly",width=13,font=self.F["ui"])
        cb.pack(side=tk.LEFT,padx=2); cb.bind("<<ComboboxSelected>>",lambda _:self._upd_slice())
        ttk.Separator(cf,orient="vertical").pack(side=tk.LEFT,fill=tk.Y,padx=6,pady=2)
        self._check(cf,"Regions",self.show_regions,self._redraw_solution,"#f0abfc")
        sf=self._strip(self.ts)
        tk.Label(sf,text="Position",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=4)
        self.sl_scale=ttk.Scale(sf,from_=0,to=30,orient=tk.HORIZONTAL,variable=self.sl_pos,
                                command=lambda _:self._debounce('_sl_after_id',self._upd_slice))
        self.sl_scale.pack(side=tk.LEFT,fill=tk.X,expand=True,padx=6)
        self.sl_lbl=tk.Label(sf,text="",bg=C["panel"],fg=C["yellow"],font=self.F["mono"],width=22,anchor="w")
        self.sl_lbl.pack(side=tk.LEFT,padx=4)
        self.fsl,self.cvsl=self._mfig(self.ts,fs=(15,7))

    def _slider_row(self,parent,var,cb,name):
        sf=self._strip(parent)
        tk.Label(sf,text=name,bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=4)
        lbl=tk.Label(sf,text="—",bg=C["panel"],fg=C["yellow"],font=self.F["mono"],width=14,anchor="w")
        lbl.pack(side=tk.LEFT,padx=2)
        sc=ttk.Scale(sf,from_=0,to=30,orient=tk.HORIZONTAL,variable=var,command=cb)
        sc.pack(side=tk.LEFT,fill=tk.X,expand=True,padx=4)
        return lbl,sc

    def _build_cut_controls(self,parent,axis_var,p1_var,p2_var,axis_callback,slider_callback,extra_widgets=None):
        cf=self._strip(parent)
        self._radio_axes(cf,axis_var,axis_callback)
        if extra_widgets:
            ttk.Separator(cf,orient="vertical").pack(side=tk.LEFT,fill=tk.Y,padx=6,pady=2)
            extra_widgets(cf)
        l1,s1=self._slider_row(parent,p1_var,slider_callback,"Position 1")
        l2,s2=self._slider_row(parent,p2_var,slider_callback,"Position 2")
        return l1,l2,s1,s2

    def _init_band_tab(self):
        def extras(cf):
            self._check(cf,"Ef (combined)",self.show_ef,self._upd_band_plot,C["yellow"])
            self._check(cf,"Efn / Efp",self.show_efnp,self._upd_band_plot,C["Efn"])
            self._check(cf,"n, p overlay",self.show_carriers,self._upd_band_plot)
            self._check(cf,"Regions",self.show_regions,self._redraw_solution,"#f0abfc")
        self.bd_lbl1,self.bd_lbl2,self.bd_scale1,self.bd_scale2=self._build_cut_controls(
            self.tb,self.bd_axis,self.bd_p1,self.bd_p2,self._upd_band_axis,
            lambda _:self._debounce('_bd_after_id',self._upd_band_plot),extras)
        self.fb,self.cvb=self._mfig(self.tb)

    def _init_curr_tab(self):
        def extras(cf):
            self._check(cf,"2-D arrows",self.cur_show_arrows2d,self._upd_curr_plot,C["green"])
            self._check(cf,"3-D arrows",self.cur_show_arrows3d,self._upd_curr_plot,C["yellow"])
            self._check(cf,"Regions",self.show_regions,self._redraw_solution,"#f0abfc")
        self.cur_lbl1,self.cur_lbl2,self.cur_scale1,self.cur_scale2=self._build_cut_controls(
            self.tc,self.cur_axis,self.cur_p1,self.cur_p2,self._upd_curr_axis,
            lambda _:self._debounce('_cur_after_id',self._upd_curr_plot),extras)
        self.fc,self.cvc=self._mfig(self.tc)

    def _init_fi_tab(self):
        def extras(cf):
            self._check(cf,"E-field",self.fi_show_field,self._upd_fi_plot,C["accent"])
            self._check(cf,"Mobility",self.fi_show_mob,self._upd_fi_plot,C["green"])
            self._check(cf,"Regions",self.show_regions,self._redraw_solution,"#f0abfc")
        self.fi_lbl1,self.fi_lbl2,self.fi_scale1,self.fi_scale2=self._build_cut_controls(
            self.tf,self.fi_axis,self.fi_p1,self.fi_p2,self._upd_fi_axis,
            lambda _:self._debounce('_fi_after_id',self._upd_fi_plot),extras)
        self.ff,self.cvf=self._mfig(self.tf)

    def _init_iv_tab(self):
        self.fiv,self.cviv=self._mfig(self.ti)

    def _init_results_tab(self):
        pw=ttk.PanedWindow(self.tr,orient=tk.VERTICAL); pw.pack(fill=tk.BOTH,expand=True)
        top=tk.Frame(pw,bg=C["panel"])
        self.res_title=tk.Label(top,text="No solution yet - press ▶ Run (F5)",bg=C["panel"],fg=C["accent"],
                                font=self.F["head"],anchor="w",justify="left")
        self.res_title.pack(fill=tk.X,padx=8,pady=(6,2))
        self.res_title.bind("<Configure>",lambda e:self.res_title.configure(wraplength=max(200,e.width-10)))
        cols=("Terminal","Net","BC","V applied","V solved","I (A)","J (A/cm²)","Resolution (A)","Area (cm²)","Note")
        tf,self.res_tree=self._tree(top,cols,[170,70,60,80,80,110,110,110,100,220],height=7,fixed=False)
        for c2 in cols[3:9]: self.res_tree.column(c2,anchor="e")
        self.res_tree.column("Terminal",anchor="w"); self.res_tree.column("Note",anchor="w")
        tf.pack(fill=tk.BOTH,expand=True,padx=6,pady=4)
        mid=tk.Frame(pw,bg=C["panel"])
        tk.Label(mid,text="Extracted parameters (last sweep)",bg=C["panel"],fg=C["accent"],
                 font=self.F["uib"],anchor="w").pack(fill=tk.X,padx=8,pady=(4,0))
        self.res_par=self._textbox(mid,6)
        bot=tk.Frame(pw,bg=C["panel"])
        hb=tk.Frame(bot,bg=C["panel"]); hb.pack(fill=tk.X)
        tk.Label(hb,text="Run log",bg=C["panel"],fg=C["accent"],font=self.F["uib"],anchor="w").pack(side=tk.LEFT,padx=8,pady=(4,0))
        self._btn(hb,"Clear",lambda:(self.res_log.configure(state="normal"),self.res_log.delete("1.0","end"),
                                     self.res_log.configure(state="disabled")),bg=C["border"]).pack(side=tk.RIGHT,padx=6,pady=2)
        self.res_log=self._textbox(bot,8)
        pw.add(top,weight=3); pw.add(mid,weight=1); pw.add(bot,weight=2)

    def _textbox(self,parent,h):
        fr=tk.Frame(parent,bg=C["panel"]); fr.pack(fill=tk.BOTH,expand=True,padx=6,pady=4)
        t=tk.Text(fr,height=h,width=40,wrap="none",bg="#07111e",fg=C["text"],relief="flat",font=self.F["monos"],
                  padx=6,pady=4,highlightthickness=0,state="disabled")
        vs=ttk.Scrollbar(fr,orient="vertical",command=t.yview); hs=ttk.Scrollbar(fr,orient="horizontal",command=t.xview)
        t.configure(yscrollcommand=vs.set,xscrollcommand=hs.set)
        t.grid(row=0,column=0,sticky="nsew"); vs.grid(row=0,column=1,sticky="ns"); hs.grid(row=1,column=0,sticky="ew")
        fr.rowconfigure(0,weight=1); fr.columnconfigure(0,weight=1)
        return t
    def _settext(self,t,s):
        t.configure(state="normal"); t.delete("1.0","end"); t.insert("end",s); t.configure(state="disabled")
    def _log(self,msg):
        try:
            self.res_log.configure(state="normal")
            self.res_log.insert("end",time.strftime("%H:%M:%S  ")+msg+"\n")
            self.res_log.see("end"); self.res_log.configure(state="disabled")
        except Exception: pass

    def _placeholder(self,fig,cv,msg):
        fig.clear()
        fig.text(0.5,0.5,msg,ha="center",va="center",color=C["sub"],fontsize=11)
        cv.draw_idle()
    def _clear_solution_plots(self):
        msg="No solution for this device yet - press ▶ Run (F5)"
        for fig,cv in [(self.fb,self.cvb),(self.fc,self.cvc),(self.ff,self.cvf),(self.fsl,self.cvsl)]:
            self._placeholder(fig,cv,msg)
        self._placeholder(self.fiv,self.cviv,"No sweep yet - set it up in the Sweep page and press ↔ Sweep (F6)")
        self.sl_lbl.config(text="")

    # ─── TEMPLATE LOADING ──────────────────────────────────────
    def _load(self,name):
        if self.busy:
            messagebox.showinfo("Busy","A solve is running - stop it first (Esc)."); return
        try: t=TEMPLATES[name]()
        except Exception:
            messagebox.showerror("Template error",traceback.format_exc()); return
        self.tname=name
        self.tmeta={k:copy.deepcopy(t.get(k)) for k in ("nodes","clear","refine","view","ivnet","ivfloat","ivout","wl","extrude")}
        self.mesh_type.set(t.get("mesh","rect")); self.prism_ax.set("auto")
        self.Lx.set(t["Lx"]); self.Ly.set(t["Ly"]); self.Lz.set(t["Lz"])
        self.Nx.set(t["Nx"]); self.Ny.set(t["Ny"]); self.Nz.set(t["Nz"])
        self.use_graded.set(True)
        self.regs=t["regs"]; self.cons=t["cons"]; self.sheets=t.get("sheets",[])
        self.dense_x=list(t.get("dense_x",[])); self.dense_y=list(t.get("dense_y",[])); self.dense_z=list(t.get("dense_z",[]))
        self.ivS.set(t["ivS"]); self.ivE.set(t["ivE"]); self.ivN.set(t["ivN"])
        self.bv_thr.set(t.get("bv_thr",0.))
        mdl=int(t.get("models",M_BGN|M_TAU|M_VSAT))
        for v,b in ((self.m_bgn,M_BGN),(self.m_tau,M_TAU),(self.m_vsat,M_VSAT),(self.m_qc,M_QC),
                    (self.m_lb,M_LB),(self.m_tun,M_TUN)): v.set(bool(mdl&b))
        groups=con_groups(self.cons)
        tgt=t.get("ivnet")
        if not tgt:
            ivc=t.get("ivc",0); tgt=next((k for k,idx in groups if ivc in idx),groups[0][0] if groups else "")
        self._refresh_groups(tgt=tgt,flt=t.get("ivfloat") or "(none)",out=t.get("ivout") or "auto")
        self._rfr(redraw=False); self._rfc(redraw=False)
        cat=next((c for c,ns in template_groups() if name in ns),"")
        self.t_name.config(text=name); self.t_cat.config(text=cat)
        self.t_doc.configure(state="normal"); self.t_doc.delete("1.0","end")
        self.t_doc.insert("end",template_doc(name)); self.t_doc.configure(state="disabled")
        self.tmpl_btn.config(text=name if len(name)<=30 else name[:29]+"…")
        self.sol=None; self.ivd=None
        view=self.tmeta.get("view")
        if view: self.sec_plane.set(view[0]); self._sec_range(view[1])
        else: self.sec_plane.set("XY"); self._sec_range(None)
        self._show_sheets(); self._clear_solution_plots(); self._fill_results(None)
        for l in (self.bd_lbl1,self.bd_lbl2,self.cur_lbl1,self.cur_lbl2,self.fi_lbl1,self.fi_lbl2): l.configure(text="—")
        self._fresh_cuts=True
        cut=t.get("cut")                  # line-cut direction: along the device's junction sequence
        if cut not in ("X","Y","Z"):
            jx,jy,_=detect_junctions(self.regs,t["Lx"],t["Ly"],t["Lz"])
            nx_=len(set(jx)|set(self.dense_x)); ny_=len(set(jy)|set(self.dense_y))
            cut="Y" if ny_>nx_ else "X"
        for v in (self.bd_axis,self.cur_axis,self.fi_axis): v.set(cut)
        self._mesh_update()
        self.sv.set(f"{name}: {t.get('desc','')}")
        self._log(f"Loaded template '{name}'")

    def _show_sheets(self):
        if not self.sheets: self.sheet_lbl.config(text="none"); return
        self.sheet_lbl.config(text="\n".join(f"{s.label}: σ = {s.sigma:+.2e} q/cm² on {s.axis} = {s.pos:g} nm"
                                             for s in self.sheets))

    def _refresh_groups(self,tgt=None,flt=None,out=None):
        """Rebuild the sweep-terminal lists from the contacts (nets merged)."""
        groups=con_groups(self.cons)
        disp=[group_label(self.cons,k,idx) for k,idx in groups]
        self._gmap={d:k for d,(k,_) in zip(disp,groups)}
        rev={k:d for d,k in self._gmap.items()}
        cur_t=self._gmap.get(self.iv_tgt.get()) if tgt is None else tgt
        cur_f=self._gmap.get(self.iv_flt.get(),"(none)") if flt is None else flt
        cur_o=self._gmap.get(self.iv_out.get(),"auto") if out is None else out
        self.cb_tgt.configure(values=disp)
        self.cb_flt.configure(values=["(none)"]+disp)
        self.cb_out.configure(values=["auto"]+disp)
        self.iv_tgt.set(rev.get(cur_t,disp[0] if disp else ""))
        self.iv_flt.set(rev.get(cur_f,"(none)"))
        self.iv_out.set(rev.get(cur_o,"auto"))

    # ─── MESH ──────────────────────────────────────────────────
    def _device_dict(self,deep=True):
        cp=copy.deepcopy if deep else (lambda x:x)
        return dict(Lx=float(self.Lx.get()),Ly=float(self.Ly.get()),Lz=float(self.Lz.get()),
                    Nx=int(self.Nx.get()),Ny=int(self.Ny.get()),Nz=int(self.Nz.get()),
                    regs=cp(self.regs),cons=cp(self.cons),sheets=cp(self.sheets),
                    dense_x=list(self.dense_x),dense_y=list(self.dense_y),dense_z=list(self.dense_z),
                    nodes=cp(self.tmeta.get("nodes")),clear=cp(self.tmeta.get("clear")),
                    refine=cp(self.tmeta.get("refine")),
                    mesh=self.mesh_type.get(),prism=None if self.prism_ax.get()=="auto" else self.prism_ax.get(),
                    extrude=self.tmeta.get("extrude"))
    def _mesh(self):
        """Node coordinates (nm) of the present device, cached."""
        try: t=self._device_dict(deep=False)
        except Exception: return None
        key=repr((t["Lx"],t["Ly"],t["Lz"],t["Nx"],t["Ny"],t["Nz"],self.use_graded.get(),
                  [(r.mat,r.x0,r.y0,r.z0,r.x1,r.y1,r.z1,r.dtype,r.doping,_shape_key(r)) for r in t["regs"]],
                  [(c.face,c.i0p,c.i1p,c.j0p,c.j1p,c.k0p,c.k1p,c.bc) for c in t["cons"]],
                  [(s.axis,s.pos) for s in t["sheets"]],t["dense_x"],t["dense_y"],t["dense_z"],
                  repr(t["nodes"]),repr(t.get("clear")),repr(t["refine"])))
        if self._mesh_cache[0]==key: return self._mesh_cache[1]
        try:
            xs,ys,zs=make_grid(t,self.use_graded.get())
            m=(np.asarray(xs)*1e9,np.asarray(ys)*1e9,np.asarray(zs)*1e9)
        except Exception: m=None
        self._mesh_cache=(key,m); return m
    def _prism(self):
        """Triangular-prism mesh of the present device (cached); None on the
        rectangular mesh or if it cannot be built."""
        if self.mesh_type.get()!="tri": return None
        try: t=self._device_dict(deep=False)
        except Exception: return None
        key=device_signature(t,(),0,0.,True)
        if self._prism_cache[0]==key: return self._prism_cache[1]
        try:
            xs,ys,zs=make_grid(t,True)
            pm=SM.PrismMesh(t,(xs,ys,zs),prism_axis(t,t.get("prism")),INSULATORS)
        except Exception:
            traceback.print_exc(); pm=None
        self._prism_cache=(key,pm); return pm
    def _node_count(self):
        """Nodes of the mesh the next run will use."""
        if self.mesh_type.get()=="tri":
            pm=self._prism(); return pm.N if pm is not None else 0
        m=self._mesh(); return len(m[0])*len(m[1])*len(m[2]) if m else 0
    def _mesh_changed(self):
        self._mesh_update(); self._log(f"Mesh: {'triangular prisms' if self.mesh_type.get()=='tri' else 'rectangular'}"
                                       +(f", prism axis {self.prism_ax.get()}" if self.mesh_type.get()=='tri' else ""))
    def _mesh_update(self,redraw=True):
        m=self._mesh()
        if m is None: self.mesh_lbl.config(text="mesh: invalid domain / node counts"); return
        nx,ny,nz=len(m[0]),len(m[1]),len(m[2]); N=nx*ny*nz; mb=N*BYTES_PER_NODE/1024**2
        mem=f"{mb/1024:.2f} GB" if mb>=1024 else f"{mb:.0f} MB"
        h=lambda a:np.min(np.diff(a)) if len(a)>1 else 0.
        pm=self._prism()
        if self.mesh_type.get()=="tri" and pm is not None:
            st=pm.stats(); mb2=st["N"]*BYTES_PER_NODE*1.3/1024**2
            mem2=f"{mb2/1024:.2f} GB" if mb2>=1024 else f"{mb2:.0f} MB"
            self.mesh_lbl.config(text=f"prisms: {st['N2']:,} nodes x {st['NL']} layers (axis {st['axis'].upper()}) = {st['N']:,} nodes   ~{mem2}\n"
                                      f"{st['ntri']:,} triangles, min angle {st['min_angle']:.1f}° "
                                      f"(rectangular: {N:,} nodes)")
        elif self.mesh_type.get()=="tri":
            self.mesh_lbl.config(text="triangular mesh could not be built - see the terminal")
        else:
            self.mesh_lbl.config(text=f"{nx} × {ny} × {nz} = {N:,} nodes   ~{mem}\n"
                                      f"min spacing  x {h(m[0]):.2g}  y {h(m[1]):.2g}  z {h(m[2]):.2g} nm")
        auto=detect_junctions(self.regs,self.Lx.get(),self.Ly.get(),self.Lz.get()) if self.use_graded.get() else ([],[],[])
        fmt=lambda a,b:", ".join(f"{v:g}" for v in sorted(set(list(a)+list(b)))) or "-"
        self.dense_lbl.config(text=f"Dense planes (nm)  x: {fmt(self.dense_x,auto[0])}   "
                                   f"y: {fmt(self.dense_y,auto[1])}   z: {fmt(self.dense_z,auto[2])}")
        if redraw: self._ddev()

    # ─── TABLES ────────────────────────────────────────────────
    def _rfr(self,redraw=True):
        sel=self.rtree.selection(); idx=self.rtree.index(sel[0]) if sel else None
        for it in self.rtree.get_children(): self.rtree.delete(it)
        for k,r in enumerate(self.regs):
            ins=r.mat in INSULATORS
            self.rtree.insert("","end",values=(k,r.label+(" ◯" if r.shape is not None else ""),r.mat,"-" if ins else r.dtype.upper(),
                "-" if ins else f"{r.doping:.2e}",f"{r.x0:g}",f"{r.x1:g}",f"{r.y0:g}",f"{r.y1:g}",f"{r.z0:g}",f"{r.z1:g}"))
        if idx is not None and idx<len(self.regs):
            self._select(self.rtree,idx,load=False)
        if redraw: self._geometry_changed()
    def _con_range(self,c):
        e=con_extent(c,float(self.Lx.get()),float(self.Ly.get()),float(self.Lz.get()))
        parts=[]; fixed=[]
        for ax in ("x","y","z"):
            lo,hi=e[ax]
            if abs(hi-lo)<1e-9: fixed.append(f"{ax}={lo:g}")
            else: parts.append(f"{ax} {lo:g}–{hi:g}")
        return ", ".join(parts)+(("  @ "+", ".join(fixed)) if fixed else "")
    def _rfc(self,redraw=True):
        sel=self.ctree.selection(); idx=self.ctree.index(sel[0]) if sel else None
        for it in self.ctree.get_children(): self.ctree.delete(it)
        for k,c in enumerate(self.cons):
            gate=c.bc==2 and c.face!=6
            tx=(f"{c.eot():.3g}"+(f" ({Con.stack_str(c.stack)})" if c.stack else "")) if gate else "-"
            self.ctree.insert("","end",values=(k,c.label+(" ◯" if getattr(c,"shape",None) is not None else ""),
                c.net or "",Con.FNAME[c.face],Con.BCNAME[c.bc],
                f"{c.V:g}",f"{c.pm:g}" if c.bc else "-",tx,
                f"{c.nox:.1e}" if gate and c.nox else "-",self._con_range(c)))
        if idx is not None and idx<len(self.cons):
            self._select(self.ctree,idx,load=False)
        self._refresh_groups()
        if redraw: self._geometry_changed()
    def _select(self,tv,idx,load=True):
        ch=tv.get_children()
        if 0<=idx<len(ch):
            self._noload=not load
            tv.selection_set(ch[idx]); tv.see(ch[idx])
            self.after_idle(lambda:setattr(self,"_noload",False))
    def _geometry_changed(self):
        self._debounce("_geo_after_id",lambda:self._mesh_update(),250)

    # ─── REGION CRUD ───────────────────────────────────────────
    def _reg_from_form(self):
        try:
            d=float(self.rb_dp.get())
            vals=[float(v.get()) for v in (self.rb_x0,self.rb_y0,self.rb_z0,self.rb_x1,self.rb_y1,self.rb_z1)]
        except Exception:
            messagebox.showerror("Region","Doping and coordinates must be numbers."); return None
        x0,y0,z0,x1,y1,z1=vals
        if not (x1>x0 and y1>y0 and z1>z0):
            messagebox.showerror("Region","Each upper coordinate must exceed the lower one."); return None
        dtype=self.rb_dt.get()
        if dtype=="i": dtype="n"; d=1e10
        lbl=self.rb_lbl.get().strip() or f"{self.rb_mat.get()} {dtype} {d:.0e}"
        return Reg(self.rb_mat.get(),x0,y0,z0,x1,y1,z1,dtype,d,lbl)
    def _sel_idx(self,tv):
        sel=tv.selection(); return tv.index(sel[0]) if sel else None
    def _addreg(self):
        r=self._reg_from_form()
        if r: self.regs.append(r); self._rfr(); self._select(self.rtree,len(self.regs)-1,False)
    def _editreg(self):
        i=self._sel_idx(self.rtree)
        if i is None: messagebox.showinfo("Region","Select a region in the table first."); return
        r=self._reg_from_form()
        if r:
            r.shape=self.regs[i].shape            # an outline is kept (template-defined)
            self.regs[i]=r; self._rfr()
    def _dupreg(self):
        i=self._sel_idx(self.rtree)
        if i is None: return
        r=copy.deepcopy(self.regs[i]); r.label+=" copy"; self.regs.insert(i+1,r); self._rfr(); self._select(self.rtree,i+1,False)
    def _fillreg(self):
        self.rb_x0.set(0.); self.rb_y0.set(0.); self.rb_z0.set(0.)
        self.rb_x1.set(self.Lx.get()); self.rb_y1.set(self.Ly.get()); self.rb_z1.set(self.Lz.get())
    def _loadreg(self):
        if getattr(self,"_noload",False): return
        i=self._sel_idx(self.rtree)
        if i is None or i>=len(self.regs): return
        r=self.regs[i]
        self.rb_mat.set(r.mat); self.rb_dt.set("i" if r.doping<=1e11 else r.dtype)
        self.rb_dp.set(f"{r.doping:.3g}"); self.rb_lbl.set(r.label)
        for v,x in zip((self.rb_x0,self.rb_y0,self.rb_z0,self.rb_x1,self.rb_y1,self.rb_z1),
                       (r.x0,r.y0,r.z0,r.x1,r.y1,r.z1)): v.set(x)
    def _delreg(self):
        i=self._sel_idx(self.rtree)
        if i is not None: del self.regs[i]; self._rfr()
    def _regup(self):
        i=self._sel_idx(self.rtree)
        if i: self.regs[i-1],self.regs[i]=self.regs[i],self.regs[i-1]; self._rfr(); self._select(self.rtree,i-1,False)
    def _regdn(self):
        i=self._sel_idx(self.rtree)
        if i is not None and i<len(self.regs)-1:
            self.regs[i],self.regs[i+1]=self.regs[i+1],self.regs[i]; self._rfr(); self._select(self.rtree,i+1,False)
    def _clrreg(self):
        if self.regs and messagebox.askyesno("Regions","Delete all regions?"): self.regs.clear(); self._rfr()

    # ─── CONTACT CRUD ──────────────────────────────────────────
    def _full_face(self):
        for v,x in ((self.cb_i0,0.),(self.cb_i1,100.),(self.cb_j0,0.),(self.cb_j1,100.),(self.cb_k0,0.),(self.cb_k1,100.)): v.set(x)
    def _upd_con_axes(self,*_):
        f={n:k for k,n in enumerate(Con.FNAME)}.get(self.cb_face.get(),0)
        a=Con.FAXES[f]
        if hasattr(self,"con_axis_lbl"):
            self.con_axis_lbl.config(text=(f"Box: i → X, j → Y, k → Z (% of Lx, Ly, Lz)" if f==6 else
                                           f"{Con.FNAME[f]} face: i → {a[0]}, j → {a[1]} (% of the face)"))
            self.ci_lbl.config(text=f"{a[0] if f!=6 else 'X'} from / to %")
            self.cj_lbl.config(text=f"{a[1] if f!=6 else 'Y'} from / to %")
            self.ck_lbl.config(text="Z from / to %" if f==6 else "(Box only) k %",
                               fg=C["sub"] if f==6 else C["border"])
    def _con_from_form(self):
        FM={n:k for k,n in enumerate(Con.FNAME)}; BM={n:k for k,n in enumerate(Con.BCNAME)}
        try:
            face=FM[self.cb_face.get()]; bc=BM[self.cb_bc.get()]
            vals=[float(v.get()) for v in (self.cb_V,self.cb_pm,self.cb_i0,self.cb_i1,self.cb_j0,self.cb_j1,
                                          self.cb_k0,self.cb_k1,self.cb_tox,self.cb_nox)]
        except Exception:
            messagebox.showerror("Contact","All contact fields must be numbers."); return None
        V,pm,i0,i1,j0,j1,k0,k1,tox,nox=vals
        if not (i1>i0 and j1>j0 and (face!=6 or k1>k0)):
            messagebox.showerror("Contact","Each 'to' percentage must exceed its 'from'."); return None
        try: stack=Con.parse_stack(self.cb_stack.get())
        except ValueError as e:
            messagebox.showerror("Contact",f"Gate stack: {e}.\nUse e.g. 'SiO2 0.5, HfO2 2.5' "
                                 f"(nm, semiconductor side first; materials: {', '.join(sorted(INSULATORS))})."); return None
        lbl=self.cb_lbl.get().strip() or f"{Con.FNAME[face]} {Con.BCNAME[bc]} {V:g}V"
        c=Con(face,i0/100.,i1/100.,j0/100.,j1/100.,bc,V,pm,lbl,tox,nox,k0/100.,k1/100.,
              net=self.cb_net.get().strip(),stack=stack)
        if stack: c.tox=round(c.eot(),4)
        return c
    def _sync_net(self,c):
        """All members of c's net take c's voltage."""
        if not c.net: return
        n=0
        for o in self.cons:
            if o is not c and o.net==c.net and o.V!=c.V: o.V=c.V; n+=1
        if n: self._log(f"Net '{c.net}': {n} other contact(s) set to {c.V:g} V")
    def _addcon(self):
        c=self._con_from_form()
        if c: self.cons.append(c); self._sync_net(c); self._rfc(); self._select(self.ctree,len(self.cons)-1,False)
    def _editcon(self):
        i=self._sel_idx(self.ctree)
        if i is None: messagebox.showinfo("Contact","Select a contact in the table first."); return
        c=self._con_from_form()
        if c:
            if c.face==6: c.shape=self.cons[i].shape
            self.cons[i]=c; self._sync_net(c); self._rfc()
    def _dupcon(self):
        i=self._sel_idx(self.ctree)
        if i is None: return
        c=copy.deepcopy(self.cons[i]); c.label+=" copy"; self.cons.insert(i+1,c); self._rfc(); self._select(self.ctree,i+1,False)
    def _loadcon(self):
        if getattr(self,"_noload",False): return
        i=self._sel_idx(self.ctree)
        if i is None or i>=len(self.cons): return
        c=self.cons[i]
        self.cb_face.set(Con.FNAME[c.face]); self.cb_bc.set(Con.BCNAME[c.bc])
        self.cb_V.set(c.V); self.cb_pm.set(c.pm); self.cb_net.set(c.net or "")
        for v,x in ((self.cb_i0,c.i0p),(self.cb_i1,c.i1p),(self.cb_j0,c.j0p),(self.cb_j1,c.j1p),
                    (self.cb_k0,c.k0p),(self.cb_k1,c.k1p)): v.set(round(x*100.,4))
        self.cb_tox.set(round(c.eot(),4)); self.cb_nox.set(c.nox); self.cb_lbl.set(c.label)
        self.cb_stack.set(Con.stack_str(c.stack))
    def _delcon(self):
        i=self._sel_idx(self.ctree)
        if i is not None: del self.cons[i]; self._rfc()
    def _clrcon(self):
        if self.cons and messagebox.askyesno("Contacts","Delete all contacts?"): self.cons.clear(); self._rfc()


    # ─── JOBS: run / sweep / stop ──────────────────────────────
    def _snapshot(self):
        """Everything a worker thread needs, read from the widgets NOW (the
        worker never touches Tk).  Raises ValueError with a readable message."""
        try:
            t=self._device_dict()
            solver=(int(self.max_p.get()),int(self.max_c.get()),int(self.max_g.get()),
                    float(self.tol_p.get()),float(self.tol_c.get()),1.,1.)
            T=float(self.T.get())
        except Exception:
            raise ValueError("Domain, mesh and solver fields must be numbers.")
        if not t["regs"]: raise ValueError("Add at least one region.")
        if not t["cons"]: raise ValueError("Add at least one contact.")
        if min(t["Nx"],t["Ny"],t["Nz"])<2: raise ValueError("Each base node count must be at least 2.")
        models=((M_BGN if self.m_bgn.get() else 0)|(M_TAU if self.m_tau.get() else 0)|
                (M_VSAT if self.m_vsat.get() else 0)|(M_QC if self.m_qc.get() else 0)|
                (M_LB if self.m_lb.get() else 0)|(M_TUN if self.m_tun.get() else 0))
        graded=bool(self.use_graded.get())
        N=self._node_count()
        if N:
            if N*BYTES_PER_NODE>0.8*_ram_bytes():
                raise ValueError(f"The mesh has {N:,} nodes (~{N*BYTES_PER_NODE/1024**3:.1f} GB) - more than "
                                 f"this computer's free memory allows. Reduce Nx/Ny/Nz.")
        return t,solver,models,T,graded

    def _forget(self):
        self._c_sig=None; self._c_ok=False; self._log("Solution forgotten - next run starts from equilibrium")

    def _set_threads(self):
        try: n=max(1,min(int(self.threads.get()),self.ncpu))
        except Exception: return
        self.threads.set(n); os.environ["OMP_NUM_THREADS"]=str(n)
        self._log(f"CPU threads: {n} (used from the next run)")
    def _nthreads(self):
        try: return max(1,min(int(self.threads.get()),self.ncpu))
        except Exception: return self.ncpu

    def _set_busy(self,busy):
        self.busy=busy
        st_on,st_off=("disabled","normal") if busy else ("normal","disabled")
        self.b_run.config(state=st_on); self.b_swp.config(state=st_on); self.b_stop.config(state=st_off)
        self.tmpl_btn.config(state="disabled" if busy else "normal")
        if not busy: self.pg["value"]=0.

    def _start(self,job,work,on_done):
        """Run work(job) in a thread; progress and results come back through
        a queue polled from the Tk main loop."""
        self.job=job; api_unstop(); self._q=queue.Queue(); self._set_busy(True)
        def W():
            try: res=work(job); self._q.put(("done",res))
            except MemoryError as e: self._q.put(("error",str(e)))
            except Exception: self._q.put(("error",traceback.format_exc()))
        threading.Thread(target=W,daemon=True).start()
        self._poll_id=self.after(120,lambda:self._poll(on_done))

    def _poll(self,on_done):
        job=self.job
        try:
            t,it,res,steps,phase,tot=progress()
            ph={1:"equilibrium",2:"gate ramp",3:"bias ramp"}.get(int(phase),"")
            el=time.time()-job.t0
            if job.kind=="sweep":
                frac=(job.k+min(max(t,0.),1.))/max(job.N,1)
                ex=getattr(job,"extra",0)
                txt=f"point {min(job.k+1,job.N)}/{job.N}"+(f" (+{ex} refined)" if ex else "")
            else:
                frac=min(max(t,0.),1.); txt="solving"
            if job.ev>1: txt+=f" · I=0 search #{job.ev}"
            if ph: txt+=f" · {ph} {100*t:.0f}% · Gummel {int(it)} · res {res:.1e}"
            if job.stop: txt="stopping…"
            self.pg["value"]=100.*frac; self.pg_txt.set(f"{txt} · {el:.0f} s")
        except Exception: pass
        try:
            while True:
                kind,data=self._q.get_nowait()
                if kind=="point": self._iv_point(data)
                elif kind=="done":
                    self._set_busy(False); self.pg_txt.set("")
                    try: on_done(data)
                    except Exception: messagebox.showerror("Display error",traceback.format_exc())
                    return
                elif kind=="error":
                    self._set_busy(False); self.pg_txt.set(""); self._c_sig=None
                    self.sv.set("Solver error - see the message"); self._log("ERROR: "+data.strip().splitlines()[-1])
                    messagebox.showerror("Solver error",data); return
        except queue.Empty: pass
        self._poll_id=self.after(150,lambda:self._poll(on_done))

    def _stop(self):
        if self.busy and self.job is not None:
            self.job.request_stop(); self.pg_txt.set("stopping…"); self._log("Stop requested")

    def _prepare(self,job,t,solver,models,T,graded,sig,reuse):
        """(worker) set the device up unless the C state can be continued.
        The OpenMP thread count is a per-thread setting, so it is applied
        here, in the worker thread that runs the solver."""
        api_threads(max(1,int(getattr(job,"nthreads",0) or self.ncpu)))
        if reuse and sig==self._c_sig and self._c_ok:
            return False
        self._c_ok=False
        xs,ys,zs=setup_device(t,solver=solver,models=models,T=T,graded=graded)
        self._c_sig=sig; self._grid=(xs,ys,zs)
        return True

    def _run(self):
        if self.busy: return
        try: t,solver,models,T,graded=self._snapshot()
        except ValueError as e: messagebox.showerror("Cannot run",str(e)); return
        sig=device_signature(t,solver,models,T,graded)
        cons=t["cons"]; vapp=[c.V for c in cons]
        fkey=self._gmap.get(self.iv_flt.get(),"(none)")
        flt=group_of(cons,fkey) if fkey!="(none)" else []
        job=Job("run",1); reuse=bool(self.reuse_sol.get()); job.nthreads=self._nthreads()
        nn=self._node_count()
        self.sv.set(f"Solving {self.tname} ({nn:,} nodes{', prisms' if t.get('mesh')=='tri' else ''})…")
        def work(job):
            fresh=self._prepare(job,t,solver,models,T,graded,sig,reuse)
            ok,Vx=run_point(cons,flt,job,vapp)
            self._c_ok=True
            s=pull()
            if s is not None:
                s["dev"]=dict(regs=t["regs"],cons=cons,sheets=t.get("sheets",[]),L=(t["Lx"],t["Ly"],t["Lz"]))
                s["models"]=models
            return dict(ok=ok,Vx=Vx,s=s,fresh=fresh,dt=time.time()-job.t0,flt=flt,cons=cons,vapp=vapp,
                        its=api_stat(3),stopped=job.stop)
        self._start(job,work,self._run_done)

    def _run_done(self,r):
        s=r["s"]; self.sol=s
        ok=r["ok"]; dt=r["dt"]
        state="STOPPED" if r["stopped"] else ("CONVERGED" if ok else "NOT CONVERGED (last converged bias kept)")
        extra=""
        if r["flt"] and r["Vx"] is not None:
            extra=f" · floating {r['cons'][r['flt'][0]].net or r['cons'][r['flt'][0]].label} = {r['Vx']:.4f} V"
        self.cv2.set(f"{'✓' if ok else '✗'} {dt:.1f} s · Gummel {r['its']}"+(" · from equilibrium" if r["fresh"] else " · continued"))
        self.sv.set(f"{self.tname}: {state} in {dt:.1f} s{extra}")
        self._log(f"Run {self.tname}: {state}, {dt:.1f} s, Gummel its {r['its']}"+extra)
        self._fill_results(s,r["cons"],r["vapp"],r["flt"],r.get("Vx"))
        if s:
            self._new_solution_views()
            self.nb.select(self.tb)

    def _new_solution_views(self):
        s=self.sol
        Nx,Ny,Nz=s["Nx"],s["Ny"],s["Nz"]
        keep=not getattr(self,"_fresh_cuts",True); self._fresh_cuts=False
        ax=self.sl_axis.get(); mx={"X":Nx,"Y":Ny,"Z":Nz}[ax]-1
        self.sl_scale.configure(to=mx)
        if not keep or self.sl_pos.get()>mx or self.sl_pos.get()<=0: self.sl_pos.set(mx//2)
        self._upd_band_axis(keep=keep); self._upd_curr_axis(keep=keep); self._upd_fi_axis(keep=keep)
        self._upd_slice()

    def _sweep(self):
        if self.busy: return
        try: t,solver,models,T,graded=self._snapshot()
        except ValueError as e: messagebox.showerror("Cannot sweep",str(e)); return
        try:
            Vs=float(self.ivS.get()); Ve=float(self.ivE.get()); Np=int(self.ivN.get())
            bv=float(self.bv_thr.get())
        except Exception:
            messagebox.showerror("Cannot sweep","Start, end, points and BV knee must be numbers."); return
        if not (1<=Np<=2000): messagebox.showerror("Cannot sweep","Points must be 1…2000."); return
        cons=t["cons"]; groups=con_groups(cons)
        tkey=self._gmap.get(self.iv_tgt.get())
        tgt=group_of(cons,tkey) if tkey else []
        if not tgt: messagebox.showerror("Cannot sweep","Choose the swept terminal in the Sweep page."); return
        fkey=self._gmap.get(self.iv_flt.get(),"(none)")
        flt=group_of(cons,fkey) if fkey!="(none)" else []
        if set(flt)&set(tgt): messagebox.showerror("Cannot sweep","The floating terminal cannot be the swept one."); return
        okey=self._gmap.get(self.iv_out.get(),"auto")
        sig=device_signature(t,solver,models,T,graded)
        vapp=[c.V for c in cons]
        job=Job("sweep",Np); job.nthreads=self._nthreads()
        self._live=[]; self._live_meta=dict(cons=cons,groups=groups,tkey=tkey,fkey=fkey if flt else None,
                                            okey=okey,bv=bv,Vs=Vs,Ve=Ve,wl=self.tmeta.get("wl"))
        self.sv.set(f"Sweeping {tkey}: {Vs:g} → {Ve:g} V, {Np} points"+(f", {fkey} floating" if flt else "")+"…")
        self._log(f"Sweep {self.tname}: {tkey} {Vs:g}→{Ve:g} V ({Np} pts)"+(f", floating {fkey}" if flt else ""))
        refine=bool(self.iv_refine.get()); reuse=bool(self.reuse_sol.get())
        if self.iv_live.get(): self.nb.select(self.ti)
        def work(job):
            fresh=self._prepare(job,t,solver,models,T,graded,sig,reuse)
            R=run_sweep(cons,tgt,Vs,Ve,Np,flt,job,on_point=lambda k,rec:self._q.put(("point",rec)),
                        vapp=vapp,refine=refine)
            self._c_ok=True
            A=np.array([api_conA(c) for c in range(len(cons))])
            s=pull()
            if s is not None:
                s["dev"]=dict(regs=t["regs"],cons=cons,sheets=t.get("sheets",[]),L=(t["Lx"],t["Ly"],t["Lz"]))
                s["models"]=models
            return dict(R=R,A=A,s=s,fresh=fresh,dt=time.time()-job.t0,stopped=job.stop,
                        cons=cons,vapp=vapp,tgt=tgt,flt=flt)
        self._start(job,work,self._sweep_done)

    def _iv_point(self,rec):
        self._live.append(rec)
        if not self.iv_live.get(): return
        now=time.time()
        if now-self._live_last>0.6:
            self._live_last=now
            try: self._plot_iv(live=True)
            except Exception: pass

    def _sweep_done(self,r):
        R=r["R"]; m=self._live_meta; cons=r["cons"]
        self.ivd=dict(R=R,A=r["A"],cons=cons,groups=con_groups(cons),tkey=m["tkey"],fkey=m["fkey"],
                      okey=m["okey"],bv=m["bv"],labels=[c.label for c in cons],bcs=[c.bc for c in cons],
                      vapp=r["vapp"],tname=self.tname,wl=m.get("wl"))
        n=len(R["V"]); nc=int(np.sum(R["conv"])) if n else 0
        state="STOPPED" if r["stopped"] else "done"
        self.cv2.set(f"sweep {nc}/{n} converged · {r['dt']:.1f} s")
        self.sv.set(f"Sweep {state}: {nc}/{n} points converged in {r['dt']:.1f} s (fields show the last point)")
        self._log(f"Sweep {state}: {nc}/{n} converged, {r['dt']:.1f} s")
        if r["s"] is not None:
            self.sol=r["s"]
            self._fill_results(r["s"],cons,r["vapp"],r["flt"],None,note="last sweep point")
            self._new_solution_views()
        self._plot_iv(); self.nb.select(self.ti)

    # ─── RESULTS TABLE ─────────────────────────────────────────
    def _fill_results(self,s,cons=None,vapp=None,flt=None,Vx=None,note=""):
        for it in self.res_tree.get_children(): self.res_tree.delete(it)
        if s is None:
            self.res_title.config(text="No solution yet - press ▶ Run (F5)"); return
        cons=cons or self.cons; flt=flt or []
        nconv=api_conv()
        self.res_title.config(text=f"{self.tname} - {'converged' if nconv else 'NOT converged'}"
                                   f"{' ('+note+')' if note else ''} · {_mesh_desc(s)} · "
                                   f"Kirchhoff sum {sum(s['I']):+.2e} A")
        for k,c in enumerate(cons):
            if k>=len(s["I"]): break
            I_=s["I"][k]; A_=s["A"][k]; fl=s["Ifl"][k]; Vs_=s["Vs"][k]
            J_=I_/A_*1e-4 if A_>0 else float("nan")
            va=vapp[k] if vapp else c.V
            if k in flt: va=float("nan")
            nt=[]
            tun = c.bc==2 and bool(s.get("models",0)&M_TUN)
            insul = c.bc==2 and not tun
            if insul: nt.append("insulated gate (no DC current)")
            elif tun: nt.append("gate tunnelling current")
            elif abs(I_)<fl: nt.append("below resolution")
            if k in flt: nt.append("floating: I = 0 solved")
            if c.bc!=2 and abs(Vs_-(va if va==va else Vs_))>1e-9: nt.append("bias not reached")
            self.res_tree.insert("","end",values=(f"{k}: {c.label}",c.net or "",Con.BCNAME[c.bc],
                "-" if va!=va else f"{va:.4g}",f"{Vs_:.5g}",
                "-" if insul else f"{I_:+.4e}","-" if insul or J_!=J_ else f"{J_:+.4e}",
                "-" if insul else f"{fl:.1e}",f"{A_*1e4:.3e}" if A_>0 else "-","; ".join(nt)))


    # ─── PLOT HELPERS ──────────────────────────────────────────
    def _debounce(self,attr_name,callback,delay_ms=120):
        prev=getattr(self,attr_name,None)
        if prev is not None:
            try: self.after_cancel(prev)
            except Exception: pass
        setattr(self,attr_name,self.after(delay_ms,lambda:(setattr(self,attr_name,None),callback())))

    def _redraw_all(self):
        self._ddev()
        if self.sol is not None:
            for f in (self._plot_band,self._plot_curr,self._plot_fi,self._upd_slice):
                try: f()
                except Exception: pass
        if self.ivd is not None:
            try: self._plot_iv()
            except Exception: pass

    @staticmethod
    def _box_faces(x0,x1,y0,y1,z0,z1):
        """6 quads of a box in plot coordinates (x, z, -y)."""
        P=lambda x,y,z:(x,z,-y)
        return [[P(x0,y0,z0),P(x1,y0,z0),P(x1,y1,z0),P(x0,y1,z0)],
                [P(x0,y0,z1),P(x1,y0,z1),P(x1,y1,z1),P(x0,y1,z1)],
                [P(x0,y0,z0),P(x0,y1,z0),P(x0,y1,z1),P(x0,y0,z1)],
                [P(x1,y0,z0),P(x1,y1,z0),P(x1,y1,z1),P(x1,y0,z1)],
                [P(x0,y0,z0),P(x1,y0,z0),P(x1,y0,z1),P(x0,y0,z1)],
                [P(x0,y1,z0),P(x1,y1,z0),P(x1,y1,z1),P(x0,y1,z1)]]
    @staticmethod
    def _aspect3(a3,Lx,Ly,Lz):
        m=max(Lx,Ly,Lz)
        try: a3.set_box_aspect((max(Lx,0.18*m),max(Lz,0.18*m),max(Ly,0.18*m)))
        except Exception: pass

    _CC=[C["yellow"],C["green"],"#ff5c7a",C["orange"],C["purple"],"#06b6d4","#f0abfc","#a3e635","#fb7185","#38bdf8"]
    @staticmethod
    def _reg_color(r):
        if r.mat in INSULATORS: return "#9aa7b8",0.25
        fc=MCOLS.get(r.mat,("#38bdf8","#f87171"))[0 if r.dtype=="n" else 1]
        al=float(np.clip(0.30+0.08*np.log10(max(r.doping,1e10)/1e14),0.22,0.85))
        return fc,al

    # ─── STRUCTURE ─────────────────────────────────────────────
    def _ddev(self):
        fig=self.fd; fig.clear()
        try: Lx,Ly,Lz=float(self.Lx.get()),float(self.Ly.get()),float(self.Lz.get())
        except Exception: self.cvd.draw_idle(); return
        L={"x":Lx,"y":Ly,"z":Lz}; L3=(Lx,Ly,Lz)
        plane=self.sec_plane.get(); pos=float(self.sec_pos.get())
        H,Vv,N={"XY":("x","y","z"),"YZ":("z","y","x"),"XZ":("x","z","y")}[plane]
        if self.sec_3d.get():
            gs=GridSpec(1,3,figure=fig,width_ratios=[1.0,1.55,0.62],wspace=0.08)
            a3=fig.add_subplot(gs[0],projection="3d"); a2=fig.add_subplot(gs[1]); aLg=fig.add_subplot(gs[2])
        else:
            gs=GridSpec(1,2,figure=fig,width_ratios=[2.6,0.62],wspace=0.05)
            a3=None; a2=fig.add_subplot(gs[0]); aLg=fig.add_subplot(gs[1])
        aLg.axis("off")
        a2.set_facecolor(C["bg"])
        ext=lambda r,ax:{"x":(r.x0,r.x1),"y":(r.y0,r.y1),"z":(r.z0,r.z1)}[ax]
        tol=1e-9*max(Lx,Ly,Lz); lab_pts=[]
        for r in self.regs:
            n0,n1=ext(r,N)
            if not (n0-tol<=pos<=n1+tol): continue
            fc,al=self._reg_color(r)
            box={a:ext(r,a) for a in "xyz"}
            polys=_section_polys(r.shape,box,H,Vv,N,pos) if r.shape is not None else None
            if polys is None:
                h0,h1=ext(r,H); v0,v1=ext(r,Vv)
                a2.add_patch(plt.Rectangle((h0,v0),h1-h0,v1-v0,facecolor=fc,edgecolor=C["border"],alpha=al,lw=0.6,zorder=1))
                if (h1-h0)>0.12*L[H] and (v1-v0)>0.07*L[Vv]:
                    lab_pts.append((r,(h0+h1)/2,(v0+v1)/2))
            else:
                for pg in polys:
                    a2.add_patch(_poly_patch(pg,facecolor=fc,edgecolor=C["border"],alpha=al,lw=0.6,zorder=1))
                if polys:
                    pg=max(polys,key=lambda q:np.ptp(q[0][:,0])*np.ptp(q[0][:,1]))[0]
                    if np.ptp(pg[:,0])>0.12*L[H] and np.ptp(pg[:,1])>0.07*L[Vv]:
                        lab_pts.append((r,float(pg[:,0].mean()),float(pg[:,1].mean())))
        # a region label is drawn only where no later region covers it
        for k,(r,hc,vc) in enumerate(lab_pts):
            P={H:np.array([hc]),Vv:np.array([vc]),N:np.array([pos])}
            if SM.region_index(self.regs,P["x"],P["y"],P["z"],L3)[0]!=self.regs.index(r): continue
            a2.text(hc,vc,r.label,ha="center",va="center",fontsize=6.5,color="#eef5ff",zorder=6,clip_on=True)
        from matplotlib.lines import Line2D
        hand=[]; hidden=0
        for k,c in enumerate(self.cons):
            e=con_extent(c,Lx,Ly,Lz); col=self._CC[k%len(self._CC)]
            n0,n1=e[N]
            if not (n0-tol<=pos<=n1+tol): hidden+=1; continue
            h0,h1=e[H]; v0,v1=e[Vv]
            if c.face==6 and getattr(c,"shape",None) is not None:
                for pg in _section_polys(c.shape,e,H,Vv,N,pos) or []:
                    a2.add_patch(_poly_patch(pg,facecolor=col,edgecolor=col,alpha=0.55,hatch="///",lw=1.,zorder=3))
            elif c.face==6:
                a2.add_patch(plt.Rectangle((h0,v0),h1-h0,v1-v0,facecolor=col,edgecolor=col,alpha=0.55,hatch="///",lw=1.,zorder=3))
            elif abs(h1-h0)<tol or abs(v1-v0)<tol:
                a2.plot([h0,h1],[v0,v1],color=col,lw=5,solid_capstyle="butt",zorder=4)
            else:
                a2.add_patch(plt.Rectangle((h0,v0),h1-h0,v1-v0,facecolor="none",edgecolor=col,lw=1.6,ls="--",zorder=4))
            bc={0:"ohmic",1:"Schottky",2:"gate"}[c.bc]
            hand.append(Line2D([0],[0],color=col,lw=4,label=f"{k}: {c.label} ({bc}, {c.V:g} V)"+(f" [{c.net}]" if c.net else "")))
        for sh in self.sheets:
            oth={"x":("y","z"),"y":("x","z"),"z":("x","y")}[sh.axis]
            rg={sh.axis:(sh.pos,sh.pos),oth[0]:(sh.a0,sh.a1),oth[1]:(sh.b0,sh.b1)}
            if N!=sh.axis and rg[N][0]-tol<=pos<=rg[N][1]+tol:
                a2.plot(rg[H],rg[Vv],color="#f0abfc",lw=2,ls="--",zorder=5)
                hand.append(Line2D([0],[0],color="#f0abfc",lw=2,ls="--",label=f"σ {sh.sigma:+.1e} q/cm²"))
        m=self._mesh(); pm=self._prism() if self.mesh_type.get()=="tri" else None
        if self.sec_mesh.get() and pm is not None:
            self._draw_prism_mesh(a2,pm,H,Vv,N,pos,L)
        elif m is not None and self.sec_mesh.get():
            mm={"x":m[0],"y":m[1],"z":m[2]}
            a2.vlines(mm[H],0,L[Vv],colors=C["sub"],lw=0.25,alpha=0.45,zorder=2)
            a2.hlines(mm[Vv],0,L[H],colors=C["sub"],lw=0.25,alpha=0.45,zorder=2)
        a2.set_xlim(-0.03*L[H],1.03*L[H])
        if Vv=="y": a2.set_ylim(1.03*Ly,-0.03*Ly)
        else: a2.set_ylim(-0.03*L[Vv],1.03*L[Vv])
        if 1/3<=L[H]/max(L[Vv],1e-12)<=3: a2.set_aspect("equal",adjustable="box")   # true shapes
        lab={"x":"x (nm)","y":"depth y (nm)","z":"z (nm)"}
        a2.set_xlabel(lab[H],fontsize=8,color=C["sub"]); a2.set_ylabel(lab[Vv],fontsize=8,color=C["sub"])
        a2.tick_params(labelsize=7,colors=C["sub"])
        if pm is not None: nn=f" · prism mesh {pm.N2:,}×{pm.NL} (axis {'XYZ'[pm.ax]})"
        else: nn=f" · mesh {len(m[0])}×{len(m[1])}×{len(m[2])}" if m is not None else ""
        a2.set_title(f"{plane} section at {N} = {pos:.1f} nm{nn}",fontsize=9,fontweight="bold",color=C["text"])
        if hand:
            if hidden: hand.append(Line2D([0],[0],color="none",label=f"(+{hidden} contact(s) off this plane)"))
            aLg.legend(handles=hand,loc="upper left",fontsize=7,borderaxespad=0.,framealpha=0.9,
                       title="Contacts",title_fontsize=8)
        self._sec_axes=(a2,H,Vv,N,pos)
        self._hover_txt=a2.text(0.01,0.01,"",transform=a2.transAxes,fontsize=7,color=C["yellow"],
                                va="bottom",ha="left",zorder=10,
                                bbox=dict(boxstyle="round,pad=0.25",facecolor=C["panel"],edgecolor=C["border"],alpha=0.85))
        self._hover_txt.set_visible(False)
        if not getattr(self,"_hover_cid",None):
            self._hover_cid=self.cvd.mpl_connect("motion_notify_event",self._struct_hover)
        if a3 is not None:
            a3.set_facecolor(C["bg"])
            for r in self.regs:
                if r.mat in INSULATORS: continue
                fc,al=self._reg_color(r)
                faces=(_prism_faces(r.shape,{a:ext(r,a) for a in "xyz"}) if r.shape is not None
                       else self._box_faces(r.x0,r.x1,r.y0,r.y1,r.z0,r.z1))
                a3.add_collection3d(Poly3DCollection(faces,alpha=min(al,0.35),facecolor=fc,edgecolor=C["border"],linewidth=0.25))
            for k,c in enumerate(self.cons):
                e=con_extent(c,Lx,Ly,Lz); col=self._CC[k%len(self._CC)]
                if c.face==6 and getattr(c,"shape",None) is not None:
                    faces=_prism_faces(c.shape,e,caps=not c.shape.hole); al=0.45
                else:
                    faces=self._box_faces(e["x"][0],e["x"][1],e["y"][0],e["y"][1],e["z"][0],e["z"][1]); al=0.85
                a3.add_collection3d(Poly3DCollection(faces,alpha=al,facecolor=col,edgecolor=col,linewidth=0.4))
            q={"z":[(0,0,pos),(Lx,0,pos),(Lx,Ly,pos),(0,Ly,pos)],
               "x":[(pos,0,0),(pos,Ly,0),(pos,Ly,Lz),(pos,0,Lz)],
               "y":[(0,pos,0),(Lx,pos,0),(Lx,pos,Lz),(0,pos,Lz)]}[N]
            a3.add_collection3d(Poly3DCollection([[(x,z,-y) for x,y,z in q]],alpha=0.12,facecolor="#ffffff",
                                                 edgecolor=C["yellow"],linewidth=0.8))
            a3.set_xlim(0,Lx); a3.set_ylim(0,Lz); a3.set_zlim(-Ly,0); self._aspect3(a3,Lx,Ly,Lz)
            sax3(a3,"x (nm)","z (nm)","depth y (nm)","3-D view")
            _depth_axis(a3)
        try: fig.tight_layout(pad=0.6)
        except Exception: pass
        self.cvd.draw_idle()

    def _draw_prism_mesh(self,a2,pm,H,Vv,N,pos,L):
        """Mesh overlay of a triangular-prism mesh in a section plane: the
        triangulation itself when the plane is a cross-section, else the
        layer planes and the triangle edges cut by the plane."""
        un,vn,an="xyz"[pm.iu],"xyz"[pm.iv],"xyz"[pm.ax]
        if N==an:
            hu=pm.P2[:,0] if H==un else pm.P2[:,1]; vv=pm.P2[:,1] if Vv==vn else pm.P2[:,0]
            a2.triplot(mtri.Triangulation(hu,vv,pm.T),color="#9fb3c8",lw=0.3,alpha=0.45,zorder=4.5)
            return
        # plane containing the prism axis: w = the other in-plane coordinate
        ci=0 if N==un else 1; wi=1-ci
        E=pm.tri.edges; p0=pm.P2[E[:,0]]; p1=pm.P2[E[:,1]]
        c0,c1=p0[:,ci]-pos,p1[:,ci]-pos
        cut=(c0*c1<=0)&(c0!=c1)
        t=c0[cut]/(c0[cut]-c1[cut]); w=p0[cut,wi]+t*(p1[cut,wi]-p0[cut,wi])
        on=np.abs(pm.P2[:,ci]-pos)<1e-9
        w=np.unique(np.round(np.concatenate([w,pm.P2[on,wi]]),9))
        kw=dict(colors="#9fb3c8",lw=0.3,alpha=0.4,zorder=4.5)
        if H==an:
            a2.vlines(pm.A,0,L[Vv],**kw); a2.hlines(w,0,L[H],**kw)
        else:
            a2.hlines(pm.A,0,L[H],**kw); a2.vlines(w,0,L[Vv],**kw)

    def _struct_hover(self,ev):
        sa=getattr(self,"_sec_axes",None); ht=getattr(self,"_hover_txt",None)
        if sa is None or ht is None: return
        a2,H,Vv,N,pos=sa
        if ev.inaxes is not a2 or ev.xdata is None:
            if ht.get_visible(): ht.set_visible(False); self.cvd.draw_idle()
            return
        P={H:ev.xdata,Vv:ev.ydata,N:pos}
        L=(float(self.Lx.get()),float(self.Ly.get()),float(self.Lz.get()))
        k=int(SM.region_index(self.regs,np.array([P["x"]]),np.array([P["y"]]),np.array([P["z"]]),L)[0])
        hit=self.regs[k] if k>=0 else None
        s=f"x {P['x']:.1f}  y {P['y']:.1f}  z {P['z']:.1f} nm"
        if hit is not None:
            s+=f"   {hit.label}: {hit.mat}"+("" if hit.mat in INSULATORS else f" {hit.dtype} {hit.doping:.2e} cm⁻³")
        ht.set_text(s); ht.set_visible(True); self.cvd.draw_idle()

    # ─── LINE CUTS ─────────────────────────────────────────────
    def _line_cut_at(self,arr,axis,p1,p2):
        s=self.sol; Nz,Ny,Nx=arr.shape
        if axis=="X":
            j=max(0,min(int(p1),Ny-1)); l=max(0,min(int(p2),Nz-1))
            return s["x"],arr[l,j,:],("y",s["y"][j]),("z",s["z"][l])
        if axis=="Y":
            i=max(0,min(int(p1),Nx-1)); l=max(0,min(int(p2),Nz-1))
            return s["y"],arr[l,:,i],("x",s["x"][i]),("z",s["z"][l])
        i=max(0,min(int(p1),Nx-1)); j=max(0,min(int(p2),Ny-1))
        return s["z"],arr[:,j,i],("x",s["x"][i]),("y",s["y"][j])

    def _setup_cut_sliders(self,axis_var,p1_var,p2_var,scale1,scale2,lbl1,lbl2,keep=False):
        s=self.sol; Nx,Ny,Nz=s["Nx"],s["Ny"],s["Nz"]
        ax=axis_var.get()
        s1,s2={"X":(Ny,Nz),"Y":(Nx,Nz),"Z":(Nx,Ny)}[ax]
        scale1.configure(to=max(0,s1-1)); scale2.configure(to=max(0,s2-1))
        if not (keep and 0<=p1_var.get()<s1): p1_var.set(self._default_cut(ax,1,s1))
        if not (keep and 0<=p2_var.get()<s2): p2_var.set(self._default_cut(ax,2,s2))
    def _default_cut(self,ax,which,n):
        """Initial transverse cut index: through the channel/junction region
        where it matters - just under the top surface for lateral devices."""
        s=self.sol
        other={"X":("y","z"),"Y":("x","z"),"Z":("x","y")}[ax][which-1]
        c=s[other]
        if other=="y":                    # horizontal cut: 2 nm below the top semiconductor surface
            semi=s["semi"]; rows=semi[semi.shape[0]//2,:,:].any(axis=1)
            js=np.where(rows)[0]
            if len(js):
                y0=c[js[0]]; return int(np.argmin(np.abs(c-(y0+2.))))
        return n//2
    def _upd_band_axis(self,keep=False):
        if not self.sol: return
        self._setup_cut_sliders(self.bd_axis,self.bd_p1,self.bd_p2,self.bd_scale1,self.bd_scale2,
                                self.bd_lbl1,self.bd_lbl2,keep)
        self._plot_band()
    def _upd_curr_axis(self,keep=False):
        if not self.sol: return
        self._setup_cut_sliders(self.cur_axis,self.cur_p1,self.cur_p2,self.cur_scale1,self.cur_scale2,
                                self.cur_lbl1,self.cur_lbl2,keep)
        self._plot_curr()
    def _upd_fi_axis(self,keep=False):
        if not self.sol: return
        self._setup_cut_sliders(self.fi_axis,self.fi_p1,self.fi_p2,self.fi_scale1,self.fi_scale2,
                                self.fi_lbl1,self.fi_lbl2,keep)
        self._plot_fi()
    def _upd_band_plot(self):
        if self.sol: self._plot_band()
    def _upd_curr_plot(self):
        if self.sol: self._plot_curr()
    def _upd_fi_plot(self):
        if self.sol: self._plot_fi()

    # ─── REGION OVERLAY ────────────────────────────────────────
    def _redraw_solution(self):
        if self.sol is None: return
        for f in (self._plot_band,self._plot_curr,self._plot_fi,self._upd_slice):
            try: f()
            except Exception: traceback.print_exc()
    def _owner(self,s):
        """Index of the region that owns every node of solution s (the last
        region containing it, half-open rule exactly as in setup_device);
        -1 where no region reaches.  Cached on the solution."""
        if s is None: return None
        if "owner" in s: return s["owner"]
        dev=s.get("dev")
        if dev: own=owner_grid(dev["regs"],s["x"]*1e-9,s["y"]*1e-9,s["z"]*1e-9)
        else: own=np.full((s["Nz"],s["Ny"],s["Nx"]),-1,dtype=np.int32)
        s["owner"]=own
        return own
    @staticmethod
    def _txt_px(ax,label,fs):
        """Approximate size (px) of a label at fontsize fs on this figure."""
        k=ax.figure.dpi/72.
        return 0.56*fs*k*len(label)+4., 1.25*fs*k
    def _ov_line(self,ax,coord,own,regs):
        """Shade and label the regions a line cut passes through; thin lines
        mark every region boundary (junctions, heterointerfaces, oxides).
        A label is drawn only where it fits (two staggered rows)."""
        if not self.show_regions.get() or own is None or len(coord)<2 or not regs: return
        coord=np.asarray(coord,float); n=len(coord)
        mids=np.concatenate(([coord[0]],0.5*(coord[1:]+coord[:-1]),[coord[-1]]))
        x0,x1=ax.get_xlim(); wpx=max(ax.bbox.width,1.); ppd=wpx/((x1-x0) or 1.)
        tr=ax.get_xaxis_transform(); k=0; last=[-1e30,-1e30]; row=0
        while k<n:
            j=k
            while j+1<n and own[j+1]==own[k]: j+=1
            idx=int(own[k])
            if 0<=idx<len(regs):
                r=regs[idx]; fc,_=self._reg_color(r); a,b=mids[k],mids[j+1]
                ax.axvspan(a,b,color=fc,alpha=0.12,lw=0,zorder=0)
                if k>0: ax.axvline(a,color=fc,lw=0.7,alpha=0.55,zorder=0.5)
                lab=r.label; tw,_=self._txt_px(ax,lab,6.)
                cpx=(0.5*(a+b)-x0)*ppd; room=(b-a)*ppd
                if room>=0.5*tw:                       # fits (may overhang its span a little)
                    for rr in (row,1-row):             # first free row (no overlap with its last label)
                        if cpx-0.5*tw>last[rr]+3:
                            ax.text(0.5*(a+b),0.985-0.06*rr,lab,transform=tr,ha="center",va="top",fontsize=6,
                                    color=fc,alpha=0.95,clip_on=True,zorder=6)
                            last[rr]=cpx+0.5*tw; row=1-rr; break
            k=j+1
    @staticmethod
    def _components(mask):
        """4-connected components of a 2-D boolean mask: list of index arrays."""
        H,W=mask.shape; lab=-np.ones(mask.shape,dtype=np.int64); comps=[]
        for j0,i0 in zip(*np.nonzero(mask)):
            if lab[j0,i0]>=0: continue
            c=len(comps); stack=[(j0,i0)]; lab[j0,i0]=c; members=[]
            while stack:
                j,i=stack.pop(); members.append((j,i))
                for jj,ii in ((j-1,i),(j+1,i),(j,i-1),(j,i+1)):
                    if 0<=jj<H and 0<=ii<W and mask[jj,ii] and lab[jj,ii]<0:
                        lab[jj,ii]=c; stack.append((jj,ii))
            comps.append(np.array(members))
        return comps
    def _ov_plane(self,ax,hc,vc,own2,regs,cons=None,axes=None,pos=None,L=None):
        """Region boundaries (the cell faces between nodes of different
        regions: exactly what the solver sees), region labels and the
        contacts cut by the plane on a 2-D map.  own2: shape (len(vc), len(hc)).
        Call after the axis limits are set."""
        if not self.show_regions.get() or own2 is None or not regs: return
        from matplotlib.collections import LineCollection
        hc=np.asarray(hc,float); vc=np.asarray(vc,float)
        hm=np.concatenate(([hc[0]],0.5*(hc[1:]+hc[:-1]),[hc[-1]]))
        vm=np.concatenate(([vc[0]],0.5*(vc[1:]+vc[:-1]),[vc[-1]]))
        segs=[]
        jj,ii=np.nonzero(own2[:,1:]!=own2[:,:-1])
        if len(jj): segs.append(np.stack([np.stack([hm[ii+1],vm[jj]],1),np.stack([hm[ii+1],vm[jj+1]],1)],1))
        jj,ii=np.nonzero(own2[1:,:]!=own2[:-1,:])
        if len(jj): segs.append(np.stack([np.stack([hm[ii],vm[jj+1]],1),np.stack([hm[ii+1],vm[jj+1]],1)],1))
        if segs:
            ax.add_collection(LineCollection(np.concatenate(segs),colors="white",linewidths=0.7,alpha=0.8,zorder=4))
        # labels: largest connected piece of each region, at its node nearest
        # to the piece's centroid; skipped when it would overlap a bigger one
        area=np.outer(np.diff(vm),np.diff(hm)); tot=area.sum() or 1.
        cands=[]
        for idx in np.unique(own2):
            if idx<0 or idx>=len(regs): continue
            comps=self._components(own2==idx)
            best=max(comps,key=lambda m:area[m[:,0],m[:,1]].sum())
            a=area[best[:,0],best[:,1]]
            if a.sum()<0.02*tot: continue
            hcen=np.average(hc[best[:,1]],weights=a); vcen=np.average(vc[best[:,0]],weights=a)
            d2=((hc[best[:,1]]-hcen)/((hc[-1]-hc[0]) or 1.))**2+((vc[best[:,0]]-vcen)/((vc[-1]-vc[0]) or 1.))**2
            jb,ib=best[int(np.argmin(d2))]
            cands.append((a.sum(),hc[ib],vc[jb],regs[idx].label))
        placed=[]
        for _,h,v,lab in sorted(cands,key=lambda c:-c[0]):
            X,Y=ax.transData.transform((h,v)); tw,th=self._txt_px(ax,lab,6.)
            box=(X-tw/2,Y-th/2,X+tw/2,Y+th/2)
            if any(not (box[2]<b[0] or box[0]>b[2] or box[3]<b[1] or box[1]>b[3]) for b in placed): continue
            placed.append(box)
            ax.text(h,v,lab,fontsize=6,color="white",ha="center",va="center",zorder=5,clip_on=True,
                    bbox=dict(boxstyle="round,pad=0.15",facecolor="#000000",alpha=0.35,edgecolor="none"))
        if cons and axes and L:
            Hn,Vn,Nn=axes; tol=1e-9*max(L)
            for k,c in enumerate(cons):
                e=con_extent(c,*L); col=self._CC[k%len(self._CC)]
                if not (e[Nn][0]-tol<=pos<=e[Nn][1]+tol): continue
                h0,h1=e[Hn]; v0,v1=e[Vn]
                if c.face==6 and getattr(c,"shape",None) is not None:
                    for pg in _section_polys(c.shape,e,Hn,Vn,Nn,pos) or []:
                        ax.add_patch(_poly_patch(pg,facecolor="none",edgecolor=col,lw=1.2,hatch="///",alpha=0.9,zorder=6))
                elif c.face==6:
                    ax.add_patch(plt.Rectangle((h0,v0),h1-h0,v1-v0,facecolor="none",edgecolor=col,lw=1.2,
                                               hatch="///",alpha=0.9,zorder=6))
                elif abs(h1-h0)<tol or abs(v1-v0)<tol:
                    ax.plot([h0,h1],[v0,v1],color=col,lw=4,solid_capstyle="butt",zorder=6,clip_on=False)

    # ─── BANDS ─────────────────────────────────────────────────
    def _plot_band(self):
        s=self.sol; fig=self.fb; fig.clear()
        if s is None: self.cvb.draw_idle(); return
        gs=GridSpec(1,2,figure=fig,wspace=0.35,width_ratios=[1.6,1.])
        ax=fig.add_subplot(gs[0]); ax3=fig.add_subplot(gs[1],projection="3d")
        ax.set_facecolor(C["bg"]); ax3.set_facecolor(C["bg"])
        A,p1,p2=self.bd_axis.get(),self.bd_p1.get(),self.bd_p2.get()
        x,Ec,(an1,av1),(an2,av2)=self._line_cut_at(s["Ec"],A,p1,p2)
        Ev=self._line_cut_at(s["Ev"],A,p1,p2)[1]; Efn=self._line_cut_at(s["Efn"],A,p1,p2)[1]
        Efp=self._line_cut_at(s["Efp"],A,p1,p2)[1]; n=self._line_cut_at(s["n"],A,p1,p2)[1]
        p=self._line_cut_at(s["p"],A,p1,p2)[1]
        self.bd_lbl1.configure(text=f"{an1} = {av1:.1f} nm"); self.bd_lbl2.configure(text=f"{an2} = {av2:.1f} nm")
        xl=f"{'depth y' if A=='Y' else A.lower()} (nm)"
        own=self._owner(s); regs=(s.get("dev") or {}).get("regs",[])
        ownl=self._line_cut_at(own,A,p1,p2)[1] if own is not None else None
        ax.plot(x,Ec,color=C["Ec"],lw=2.2,label="Ec"); ax.plot(x,Ev,color=C["Ev"],lw=2.2,label="Ev")
        ax.fill_between(x,Ev,Ec,color=C["accent"],alpha=0.05)
        if self.show_efnp.get():
            ax.plot(x,Efn,color=C["Efn"],lw=1.6,ls="--",label="Efn"); ax.plot(x,Efp,color=C["Efp"],lw=1.6,ls="--",label="Efp")
        if self.show_ef.get():
            wn=n/(n+p+1e-30); ax.plot(x,wn*Efn+(1.-wn)*Efp,color=C["yellow"],lw=2.0,ls="-.",label="Ef (combined)")
        if s.get("models",0)&M_QC and "Lqn" in s:
            Lqn=self._line_cut_at(s["Lqn"],A,p1,p2)[1]; Lqp=self._line_cut_at(s["Lqp"],A,p1,p2)[1]
            if np.nanmax(np.abs(np.nan_to_num(Lqn)))>1e-4 or np.nanmax(np.abs(np.nan_to_num(Lqp)))>1e-4:
                ax.plot(x,Ec+Lqn,color=C["Ec"],lw=1.1,ls=(0,(1,1)),label="Ec + Λn (quantum)")
                ax.plot(x,Ev-Lqp,color=C["Ev"],lw=1.1,ls=(0,(1,1)),label="Ev − Λp (quantum)")
        ax.legend(loc="best",fontsize=7,framealpha=0.85)
        sax(ax,xl,"Energy (eV)",f"Band diagram along {A}  ({an1} = {av1:.1f}, {an2} = {av2:.1f} nm)")
        self._ov_line(ax,x,ownl,regs)
        if self.show_carriers.get():
            a2=ax.twinx()
            a2.semilogy(x,n,color=C["Ec"],lw=1.2,ls=":",alpha=0.85,label="n")
            a2.semilogy(x,p,color=C["Ev"],lw=1.2,ls=":",alpha=0.85,label="p")
            a2.set_ylabel("n, p (cm⁻³)",fontsize=8,color=C["sub"]); a2.tick_params(labelsize=7,colors=C["sub"])
        for d in {"X":self.dense_x,"Y":self.dense_y,"Z":self.dense_z}[A]:
            ax.axvline(d,color=C["yellow"],lw=0.8,ls=":",alpha=0.35)
        Y=np.linspace(0,0.5,5); X2,Y2=np.meshgrid(x,Y)
        ax3.plot_surface(X2,Y2,np.tile(Ec,(5,1)),color=C["Ec"],alpha=0.7,linewidth=0)
        ax3.plot_surface(X2,Y2,np.tile(Ev,(5,1)),color=C["Ev"],alpha=0.7,linewidth=0)
        if self.show_efnp.get():
            ax3.plot(x,np.full_like(x,.25),Efn,color=C["Efn"],lw=2,ls="--")
            ax3.plot(x,np.full_like(x,.25),Efp,color=C["Efp"],lw=2,ls="--")
        sax3(ax3,xl,"","E (eV)",f"Band ribbon along {A}"); ax3.set_yticks([])
        try: fig.tight_layout(pad=1.)
        except Exception: pass
        self.cvb.draw_idle()

    # ─── CURRENT ───────────────────────────────────────────────
    def _plot_curr(self):
        s=self.sol; fig=self.fc; fig.clear()
        if s is None: self.cvc.draw_idle(); return
        gs=GridSpec(2,2,figure=fig,hspace=0.55,wspace=0.42)
        aJ=fig.add_subplot(gs[0,0]); aR=fig.add_subplot(gs[1,0])
        aQ=fig.add_subplot(gs[0,1]); aQ3=fig.add_subplot(gs[1,1],projection="3d")
        for a in (aJ,aR,aQ,aQ3): a.set_facecolor(C["bg"])
        A,p1,p2=self.cur_axis.get(),self.cur_p1.get(),self.cur_p2.get()
        cx,Jx,(an1,av1),(an2,av2)=self._line_cut_at(s["Jx"],A,p1,p2)
        Jy=self._line_cut_at(s["Jy"],A,p1,p2)[1]; Jz=self._line_cut_at(s["Jz"],A,p1,p2)[1]
        Rt=self._line_cut_at(s["R"],A,p1,p2)[1]
        self.cur_lbl1.configure(text=f"{an1} = {av1:.1f} nm"); self.cur_lbl2.configure(text=f"{an2} = {av2:.1f} nm")
        xl=f"{'depth y' if A=='Y' else A.lower()} (nm)"
        own=self._owner(s); dev=s.get("dev") or {}; regs=dev.get("regs",[])
        ownl=self._line_cut_at(own,A,p1,p2)[1] if own is not None else None
        aJ.plot(cx,Jx,color="#34d399",lw=1.6,label="Jx"); aJ.plot(cx,Jy,color="#fbbf24",lw=1.6,label="Jy")
        aJ.plot(cx,Jz,color="#f87171",lw=1.6,label="Jz")
        aJ.plot(cx,np.sqrt(Jx**2+Jy**2+Jz**2),color=C["green"],lw=2.2,ls="--",label="|J|")
        aJ.axhline(0,color=C["border"],lw=0.5); aJ.legend(fontsize=7,framealpha=0.85)
        sax(aJ,xl,"J (A/cm²)",f"Current density along {A} ({an1} = {av1:.1f}, {an2} = {av2:.1f} nm)")
        if np.any(np.abs(np.nan_to_num(Rt))>1e-30):
            aR.plot(cx,Rt,color=C["red"],lw=2.,label="R net"); aR.legend(fontsize=7)
        else:
            aR.text(0.5,0.5,"R ≈ 0 (near equilibrium)",transform=aR.transAxes,ha="center",va="center",fontsize=10,color=C["sub"])
        sax(aR,xl,"R (cm⁻³s⁻¹)",f"Net recombination along {A}")
        self._ov_line(aJ,cx,ownl,regs); self._ov_line(aR,cx,ownl,regs)
        Nz_,Ny_,Nx_=s["Jx"].shape
        # in-plane map: the plane that CONTAINS the cut line
        if A in ("X","Y"):
            l=max(0,min(int(p2),Nz_-1))
            U=s["Jx"][l,:,:]; Vq=s["Jy"][l,:,:]; hx,vy=s["x"],s["y"]; xlab,ylab="x (nm)","depth y (nm)"
            ttl=f"XY plane @ z = {s['z'][l]:.1f} nm"; invert=True
            cut=("h",s["y"][max(0,min(int(p1),Ny_-1))]) if A=="X" else ("v",s["x"][max(0,min(int(p1),Nx_-1))])
            own2=own[l,:,:] if own is not None else None; ovax=("x","y","z"); ovpos=s["z"][l]
        else:
            i=max(0,min(int(p1),Nx_-1))
            U=s["Jz"][:,:,i].T; Vq=s["Jy"][:,:,i].T; hx,vy=s["z"],s["y"]; xlab,ylab="z (nm)","depth y (nm)"
            ttl=f"YZ plane @ x = {s['x'][i]:.1f} nm"; invert=True
            cut=("h",s["y"][max(0,min(int(p2),Ny_-1))])
            own2=own[:,:,i].T if own is not None else None; ovax=("z","y","x"); ovpos=s["x"][i]
        M=np.sqrt(U**2+Vq**2)
        if np.nanmax(M)>1e-30:
            X1,Y1=np.meshgrid(hx,vy)
            pc=aQ.pcolormesh(X1,Y1,np.ma.masked_invalid(M),cmap="viridis",shading="nearest")
            fig.colorbar(pc,ax=aQ,label="|J| (A/cm²)",shrink=0.85,pad=0.02)
            if self.cur_show_arrows2d.get():
                gx=np.linspace(hx[0],hx[-1],22); gy=np.linspace(vy[0],vy[-1],14)
                ii=np.clip(np.searchsorted(hx,gx),0,len(hx)-1); jj=np.clip(np.searchsorted(vy,gy),0,len(vy)-1)
                Us=U[np.ix_(jj,ii)]; Vs=Vq[np.ix_(jj,ii)]; ms=np.sqrt(Us**2+Vs**2)
                okm=ms>1e-30*np.nanmax(M)+1e-300
                ln=np.log10(np.where(okm,ms,1e-300)); lo,hi=np.nanmin(ln[okm]) if okm.any() else 0,np.nanmax(ln)
                sc=0.3+0.7*(ln-lo)/(hi-lo+1e-30); sc=np.where(okm,sc,0.)
                Gx,Gy=np.meshgrid(hx[ii],vy[jj])
                aQ.quiver(Gx,Gy,np.where(okm,Us/np.where(okm,ms,1),0)*sc,
                          np.where(okm,Vs/np.where(okm,ms,1),0)*sc*(-1 if invert else 1),
                          color="white",scale=25,width=0.004,headwidth=4,headlength=5,alpha=0.9,pivot="mid")
        else:
            aQ.text(0.5,0.5,"J ≈ 0 (no current flow)",transform=aQ.transAxes,ha="center",va="center",fontsize=10,color=C["sub"])
        if cut[0]=="h": aQ.axhline(cut[1],color=C["yellow"],lw=0.9,ls="--",alpha=0.8)
        else: aQ.axvline(cut[1],color=C["yellow"],lw=0.9,ls="--",alpha=0.8)
        aQ.set_xlim(hx[0],hx[-1]); aQ.set_ylim(vy[0],vy[-1])
        self._ov_plane(aQ,hx,vy,own2,regs,dev.get("cons"),ovax,ovpos,dev.get("L"))
        sax(aQ,xlab,ylab,ttl+"  (dashed: the line cut)")
        if invert: aQ.invert_yaxis()
        if self.cur_show_arrows3d.get():
            x,y,z=s["x"],s["y"],s["z"]
            sxk=max(1,Nx_//8); syk=max(1,Ny_//6); szk=max(1,Nz_//5)
            Xa,Ya,Za=np.meshgrid(x[::sxk],y[::syk],z[::szk],indexing="xy")
            Ua=s["Jx"][::szk,::syk,::sxk].transpose(1,2,0); Va=s["Jy"][::szk,::syk,::sxk].transpose(1,2,0)
            Wa=s["Jz"][::szk,::syk,::sxk].transpose(1,2,0)
            Ma=np.nan_to_num(np.sqrt(Ua**2+Va**2+Wa**2))
            if Ma.max()>1e-30:
                msf=np.where(Ma>1e-30,Ma,1.); Lr=0.08*max(x.max(),y.max(),z.max())
                from matplotlib.colors import LogNorm
                mf=Ma.flatten(); mpos=mf[mf>1e-30]
                norm=LogNorm(vmin=max(mpos.min(),Ma.max()*1e-4),vmax=Ma.max())
                cols=plt.cm.viridis(norm(np.clip(mf,norm.vmin,None)))
                aQ3.quiver(Xa.flatten(),Za.flatten(),-Ya.flatten(),(Ua/msf*Lr).flatten(),(Wa/msf*Lr).flatten(),
                           (-Va/msf*Lr).flatten(),colors=cols,arrow_length_ratio=0.4,linewidth=1.2)
            aQ3.set_xlim(0,x.max()); aQ3.set_ylim(0,z.max()); aQ3.set_zlim(-y.max(),0)
            self._aspect3(aQ3,x.max(),y.max(),z.max())
            sax3(aQ3,"x","z","depth y","3-D current field"); _depth_axis(aQ3)
        else:
            aQ3.text2D(0.5,0.5,"3-D arrows hidden",transform=aQ3.transAxes,ha="center",fontsize=10,color=C["sub"])
            sax3(aQ3,"x","z","depth y","3-D current field (off)")
        try: fig.tight_layout(pad=1.)
        except Exception: pass
        self.cvc.draw_idle()

    # ─── FIELD & MOBILITY ──────────────────────────────────────
    def _plot_fi(self):
        s=self.sol; fig=self.ff; fig.clear()
        if s is None: self.cvf.draw_idle(); return
        gs=GridSpec(1,2,figure=fig,wspace=0.38)
        aE=fig.add_subplot(gs[0]); aM=fig.add_subplot(gs[1])
        aE.set_facecolor(C["bg"]); aM.set_facecolor(C["bg"])
        A,p1,p2=self.fi_axis.get(),self.fi_p1.get(),self.fi_p2.get()
        cx,phi,(an1,av1),(an2,av2)=self._line_cut_at(s["phi"],A,p1,p2)
        mn=self._line_cut_at(s["mn"],A,p1,p2)[1]; mp=self._line_cut_at(s["mp"],A,p1,p2)[1]
        self.fi_lbl1.configure(text=f"{an1} = {av1:.1f} nm"); self.fi_lbl2.configure(text=f"{an2} = {av2:.1f} nm")
        xl=f"{'depth y' if A=='Y' else A.lower()} (nm)"
        own=self._owner(s); regs=(s.get("dev") or {}).get("regs",[])
        ownl=self._line_cut_at(own,A,p1,p2)[1] if own is not None else None
        E=-np.gradient(phi)/(np.gradient(cx)*1e-7)
        if self.fi_show_field.get():
            aE.plot(cx,E,color=C["accent"],lw=2,label=f"E{A.lower()}")
            aE.axhline(0,color=C["border"],lw=0.5)
            aE.fill_between(cx,E,0,where=E>0,color=C["accent"],alpha=0.15)
            aE.fill_between(cx,E,0,where=E<0,color=C["red"],alpha=0.15); aE.legend(fontsize=8)
        else:
            aE.text(0.5,0.5,"E-field hidden",transform=aE.transAxes,ha="center",va="center",fontsize=10,color=C["sub"])
        sax(aE,xl,"E (V/cm)",f"Electric field along {A} ({an1} = {av1:.1f}, {an2} = {av2:.1f} nm)")
        if self.fi_show_mob.get():
            aM.plot(cx,mn,color=C["Ec"],lw=2,label="µn"); aM.plot(cx,mp,color=C["Ev"],lw=2,label="µp"); aM.legend(fontsize=8)
        else:
            aM.text(0.5,0.5,"Mobility hidden",transform=aM.transAxes,ha="center",va="center",fontsize=10,color=C["sub"])
        sax(aM,xl,"µ (cm²/Vs)",f"Carrier mobility along {A}")
        self._ov_line(aE,cx,ownl,regs); self._ov_line(aM,cx,ownl,regs)
        try: fig.tight_layout(pad=1.)
        except Exception: pass
        self.cvf.draw_idle()

    # ─── 3-D SLICER ────────────────────────────────────────────
    def _upd_slice(self,*_):
        s=self.sol
        if s is None: return
        axis=self.sl_axis.get(); qty=self.sl_qty.get()
        Nx,Ny,Nz=s["Nx"],s["Ny"],s["Nz"]
        mx={"X":Nx-1,"Y":Ny-1,"Z":Nz-1}[axis]; self.sl_scale.configure(to=mx)
        pos=max(0,min(int(round(self.sl_pos.get())),mx))
        def qcalc(s,semi):
            if qty=="phi": return s["phi"],"RdBu_r","φ (V)"
            if qty=="log10(n)": return np.log10(s["n"]),"Blues","log₁₀ n (cm⁻³)"
            if qty=="log10(p)": return np.log10(s["p"]),"Reds","log₁₀ p (cm⁻³)"
            if qty=="Ec": return np.where(semi,s["Ec"],np.nan),"plasma","Ec (eV)"
            if qty=="Ev": return np.where(semi,s["Ev"],np.nan),"viridis","Ev (eV)"
            if qty=="Efn": return s["Efn"],"Greens","Efn (eV)"
            if qty=="Efp": return s["Efp"],"Oranges","Efp (eV)"
            if qty=="log10|J|":
                return (np.where(semi,np.log10(np.sqrt(s["Jx"]**2+s["Jy"]**2+s["Jz"]**2)+1e-30),np.nan),
                        "hot","log₁₀|J| (A/cm²)")
            if qty=="mobility µn": return s["mn"],"YlOrRd","µn (cm²/Vs)"
            if qty=="net doping":
                Nn=s["Nnet"]; return np.where(semi,np.sign(Nn)*np.log10(1.+np.abs(Nn)),np.nan),"RdBu_r","±log₁₀|Nd−Na| (n > 0)"
            if qty=="log10|J gate|":
                Jg=np.abs(s.get("Jg",np.zeros_like(s["phi"])))
                return (np.where(semi&(Jg>0),np.log10(Jg+1e-300),np.nan),"inferno",
                        "log₁₀|J gate tunnelling| (A/cm²)"+("" if self.sol.get("models",0)&M_TUN else " - tunnelling off"))
            if qty=="Λn quantum":
                return (np.where(semi,s.get("Lqn",np.zeros_like(s["phi"])),np.nan),"magma",
                        "Λn quantum potential (eV)"+("" if self.sol.get("models",0)&M_QC else " - MLDA off"))
            Rt=s["R"]
            return np.where(semi,np.where(np.abs(Rt)>1e-30,np.log10(np.abs(Rt)+1e-30),np.nan),np.nan),"PuOr","log₁₀|R| (cm⁻³s⁻¹)"
        Q3,cmap,lbl=qcalc(s,s["semi"])
        pm=s.get("prism") if s.get("mesh")=="tri" else None
        Qr=qcalc(s["raw"],s["raw"]["semi"])[0] if (pm is not None and "raw" in s) else None
        x,y,z=s["x"],s["y"],s["z"]
        sx,sy,sz=Nx//2,Ny//2,Nz//2
        if axis=="X": sx=pos
        elif axis=="Y": sy=pos
        else: sz=pos
        cval={"X":x[sx],"Y":y[sy],"Z":z[sz]}[axis]
        self.sl_lbl.config(text=f"{axis.lower()} = {cval:.1f} nm  ({pos}/{mx})")
        Zxy=Q3[sz,:,:]; Zxz=Q3[:,sy,:]; Zyz=Q3[:,:,sx]
        fin=Q3[np.isfinite(Q3)]
        vmin,vmax=(float(fin.min()),float(fin.max())) if fin.size else (0.,1.)
        if qty=="log10|J gate|": vmin=max(vmin,vmax-8.)          # the top 8 decades
        if vmin==vmax: vmax=vmin+1e-10
        norm=plt.Normalize(vmin,vmax); sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap)
        fig=self.fsl; fig.clear()
        gs=GridSpec(3,2,figure=fig,width_ratios=[1.35,1.0],hspace=0.55,wspace=0.18)
        a3=fig.add_subplot(gs[:,0],projection="3d"); a3.set_facecolor(C["bg"])
        axy=fig.add_subplot(gs[0,1]); ayz=fig.add_subplot(gs[1,1]); axz=fig.add_subplot(gs[2,1])
        for a in (axy,ayz,axz): a.set_facecolor(C["bg"])
        ds=max(1,int((Nx*Ny*Nz/40_000)**0.5)) if Nx*Ny*Nz>200_000 else 1
        xd,yd,zd=x[::ds],y[::ds],z[::ds]
        def rgba(Z):
            c=sm.to_rgba(np.nan_to_num(Z,nan=vmin)); c[...,3]=np.where(np.isfinite(Z),0.9,0.08); return c
        Xg,Yg=np.meshgrid(xd,yd)
        a3.plot_surface(Xg,np.full_like(Xg,z[sz]),-Yg,facecolors=rgba(Zxy[::ds,::ds]),shade=False,
                        rstride=max(1,len(yd)//40),cstride=max(1,len(xd)//40))
        Xg2,Zg2=np.meshgrid(xd,zd)
        a3.plot_surface(Xg2,Zg2,np.full_like(Xg2,-y[sy]),facecolors=rgba(Zxz[::ds,::ds]),shade=False,
                        rstride=max(1,len(zd)//40),cstride=max(1,len(xd)//40))
        Yg3,Zg3=np.meshgrid(yd,zd)
        a3.plot_surface(np.full_like(Yg3,x[sx]),Zg3,-Yg3,facecolors=rgba(Zyz[::ds,::ds]),shade=False,
                        rstride=max(1,len(zd)//40),cstride=max(1,len(yd)//40))
        a3.set_xlim(x[0],x[-1]); a3.set_ylim(z[0],z[-1]); a3.set_zlim(-y[-1],-y[0])
        self._aspect3(a3,x[-1]-x[0],y[-1]-y[0],z[-1]-z[0])
        sax3(a3,"x (nm)","z (nm)","depth y (nm)",f"{lbl}"); _depth_axis(a3)
        fig.colorbar(sm,ax=a3,shrink=0.55,pad=0.06,label=lbl)
        own=self._owner(s); dev=s.get("dev") or {}; regs=dev.get("regs",[])
        def pmap(ax,hc,vc,Z,title,xl,yl,hl=None,vl=None,inv=False,o2=None,ovax=None,ovpos=None):
            if Qr is not None and ovax is not None and ovax[2]=="xyz"[pm.ax]:
                # cross-section of a prism mesh: the node values on the triangulation itself
                l=int(np.argmin(np.abs(pm.A-ovpos)))
                vals=Qr[pm.N2*l:pm.N2*(l+1)]
                un="xyz"[pm.iu]
                th=pm.P2[:,0] if ovax[0]==un else pm.P2[:,1]; tv=pm.P2[:,1] if ovax[0]==un else pm.P2[:,0]
                tr=mtri.Triangulation(th,tv,pm.T)
                bad=~np.isfinite(vals)
                tr.set_mask(np.any(bad[pm.T],axis=1))
                ax.tripcolor(tr,np.where(bad,vmin,vals),cmap=cmap,vmin=vmin,vmax=vmax,shading="gouraud")
                ax.set_xlim(hc[0],hc[-1]); ax.set_ylim(vc[0],vc[-1])
                wh=(vc[-1]-vc[0])/max(hc[-1]-hc[0],1e-12)
                if 0.5<=wh<=2.:          # keep a round wire round (true aspect when the section is compact)
                    ax.set_aspect("equal",adjustable="box")
                if o2 is not None and dev:
                    hf=np.linspace(hc[0],hc[-1],400); vf=np.linspace(vc[0],vc[-1],400)
                    H2,V2=np.meshgrid(hf,vf)
                    P={ovax[0]:H2,ovax[1]:V2,ovax[2]:np.full_like(H2,ovpos)}
                    o2f=SM.region_index(regs,P["x"],P["y"],P["z"],dev.get("L"))
                    self._ov_plane(ax,hf,vf,o2f,regs,dev.get("cons"),ovax,ovpos,dev.get("L"))
            else:
                pc=ax.pcolormesh(hc,vc,np.ma.masked_invalid(Z),cmap=cmap,vmin=vmin,vmax=vmax,shading="nearest")
            if o2 is not None and not (Qr is not None and ovax is not None and ovax[2]=="xyz"[pm.ax]):
                ax.set_xlim(hc[0],hc[-1]); ax.set_ylim(vc[0],vc[-1])
                self._ov_plane(ax,hc,vc,o2,regs,dev.get("cons"),ovax,ovpos,dev.get("L"))
            if hl is not None: ax.axvline(hl,color=C["yellow"],lw=0.9,ls="--",alpha=0.8)
            if vl is not None: ax.axhline(vl,color=C["green"],lw=0.9,ls="--",alpha=0.8)
            ax.set_title(title,fontsize=8,color=C["text"])
            ax.set_xlabel(xl,fontsize=7,color=C["sub"]); ax.set_ylabel(yl,fontsize=7,color=C["sub"])
            ax.tick_params(labelsize=6,colors=C["sub"])
            if inv: ax.invert_yaxis()
        pmap(axy,x,y,Zxy,f"XY @ z = {z[sz]:.1f} nm","x (nm)","depth y (nm)",x[sx],y[sy],True,
             own[sz,:,:] if own is not None else None,("x","y","z"),z[sz])
        pmap(ayz,z,y,Zyz.T,f"YZ @ x = {x[sx]:.1f} nm","z (nm)","depth y (nm)",z[sz],y[sy],True,
             own[:,:,sx].T if own is not None else None,("z","y","x"),x[sx])
        pmap(axz,x,z,Zxz,f"XZ @ y = {y[sy]:.1f} nm","x (nm)","z (nm)",x[sx],z[sz],False,
             own[:,sy,:] if own is not None else None,("x","z","y"),y[sy])
        self.cvsl.draw_idle()


    # ─── I-V / TRANSFER CURVES ─────────────────────────────────
    def _iv_data(self,live):
        if live:
            recs=self._live; m=getattr(self,"_live_meta",None)
            if not recs or m is None: return None
            V=np.array([r["V"] for r in recs]); o=np.argsort(V,kind="stable")
            if m["Ve"]<m["Vs"]: o=o[::-1]
            R=dict(V=V[o],I=np.array([r["I"] for r in recs])[o],F=np.array([r["F"] for r in recs])[o],
                   conv=np.array([r["conv"] for r in recs],bool)[o],
                   Vf=np.array([np.nan if r["Vf"] is None else r["Vf"] for r in recs],float)[o])
            return dict(R=R,A=None,cons=m["cons"],groups=m["groups"],tkey=m["tkey"],fkey=m["fkey"],
                        okey=m["okey"],bv=m["bv"],wl=m.get("wl"))
        return self.ivd

    def _plot_iv(self,live=False):
        d=self._iv_data(live); fig=self.fiv
        if d is None or len(d["R"]["V"])==0:
            if not live: self._placeholder(fig,self.cviv,"No sweep yet - set it up in the Sweep page and press ↔ Sweep (F6)")
            return
        fig.clear()
        R=d["R"]; V=R["V"]; I=R["I"]; F=R["F"]; cv=R["conv"]
        cons=d["cons"]; groups=d["groups"]
        gI={k:I[:,idx].sum(axis=1) for k,idx in groups}
        gF={k:F[:,idx].sum(axis=1) for k,idx in groups}
        gate={k:all(cons[c].bc==2 for c in idx) for k,idx in groups}
        A=d.get("A")
        gA={k:(float(sum(A[c] for c in idx)) if A is not None else 0.) for k,idx in groups}
        tkey,fkey=d["tkey"],d["fkey"]
        cand=[k for k,_ in groups if not gate[k] and k!=fkey]
        big=lambda k:float(np.nanmax(np.abs(gI[k]))) if len(gI[k]) else 0.
        okey=d.get("okey")
        if okey in (None,"auto") or okey not in gI:
            if fkey: okey=max([k for k in cand if k!=tkey] or cand or [tkey],key=big)
            elif tkey in gI and not gate[tkey]: okey=tkey
            else: okey=max(cand or [tkey],key=big)
        par=[]
        if fkey:
            self._plot_vtc(fig,V,R["Vf"],cv,gI,gF,gate,tkey,fkey,okey,par)
        else:
            self._plot_ivc(fig,V,cv,gI,gF,gA,gate,tkey,okey,d.get("bv",0.),par,cand,d.get("wl"))
        try: fig.tight_layout(pad=1.)
        except Exception: pass
        self.cviv.draw_idle()
        if not live:
            hdr=f"{d.get('tname',self.tname)}: sweep of {tkey}"+(f", {fkey} floating" if fkey else "")
            self._settext(self.res_par,hdr+"\n"+("\n".join(par) if par else "(no parameters apply to this sweep)"))

    def _plot_ivc(self,fig,V,cv,gI,gF,gA,gate,tkey,okey,bv,par,cand,wl=None):
        gs=GridSpec(1,2,figure=fig,wspace=0.45)
        aL=fig.add_subplot(gs[0]); aS=fig.add_subplot(gs[1])
        aL.set_facecolor(C["bg"]); aS.set_facecolor(C["bg"])
        Io=gI[okey]; Fo=gF[okey]
        good=cv&(np.abs(Io)>=Fo); low=cv&~(np.abs(Io)>=Fo)
        if good.any(): aL.plot(V[good],Io[good],color=C["green"],lw=2,marker="o",ms=3.5,label="converged")
        if low.any(): aL.plot(V[low],Io[low],color=C["sub"],marker="o",mfc="none",ms=4,ls="none",label="below resolution")
        if (~cv).any(): aL.plot(V[~cv],Io[~cv],color=C["red"],marker="x",ms=6,ls="none",label="not converged")
        aL.axhline(0,color=C["border"],lw=0.5); aL.axvline(0,color=C["border"],lw=0.5)
        aL.legend(fontsize=7); sax(aL,f"V({tkey}) (V)",f"I({okey}) (A per cell)",f"{okey} current vs {tkey} voltage")
        Aok=gA.get(okey,0.)
        if Aok>0:
            a2=aL.twinx(); y0,y1=aL.get_ylim(); a2.set_ylim(y0/(Aok*1e4),y1/(Aok*1e4))
            a2.set_ylabel("J (A/cm²)",fontsize=8,color=C["sub"]); a2.tick_params(labelsize=7,colors=C["sub"])
        tun=[k for k in gI if gate.get(k) and len(gI[k]) and np.nanmax(np.abs(gI[k]))>0]
        for n,k in enumerate(cand+[k for k in tun if k not in cand]):
            Ic=np.abs(gI[k]); ok=cv&(Ic>=gF[k])&(Ic>0)
            if ok.any():
                aS.semilogy(V[ok],Ic[ok],color=TCOLS[n%len(TCOLS)],lw=2.2 if k==okey else 1.2,
                            ls="--" if k in tun else "-",marker="o" if k==okey else None,ms=3,
                            label=f"|I| {k}"+(" (tunnelling)" if k in tun else ""))
        if np.any(gF[okey]>0):
            aS.semilogy(V,np.maximum(gF[okey],1e-300),color=C["sub"],lw=0.8,ls=":",label=f"resolution ({okey})")
        aS.legend(fontsize=6); sax(aS,f"V({tkey}) (V)","|I| (A per cell)","All terminals (semi-log)")
        # ── parameters ──
        a=np.abs(Io); okp=good&(a>0)
        if gate.get(okey,False) and okp.any():
            j=int(np.argmax(a[okp])); Aa=gA.get(okey,0.)
            Jg=a[okp][j]/(Aa*1e4) if Aa>0 else float("nan")
            par.append(f"Gate tunnelling current |I_G| = {a[okp][j]:.3e} A per cell at V = {V[okp][j]:.3f} V"
                       +(f"  (J_G = {Jg:.3e} A/cm²)" if Jg==Jg else ""))
            if tkey==okey and okp.sum()>=3 and Aa>0:
                o=np.argsort(V[okp]); Vk=V[okp][o]; Jk=a[okp][o]/(Aa*1e4)
                for Vq in (0.5,1.0,-0.5,-1.0):
                    if Vk.min()<=Vq<=Vk.max():
                        par.append(f"J_G({Vq:+.1f} V) = {10**np.interp(Vq,Vk,np.log10(Jk)):.3e} A/cm²")
        if gate.get(tkey,False) and not gate.get(okey,False) and okp.sum()>=3:
            Vg=V[okp]; Ia=a[okp]
            lg=np.log10(Ia); dv=np.diff(Vg); dl=np.diff(lg)
            sub=(Ia[1:]<0.01*Ia.max())&(Ia[:-1]<0.01*Ia.max())&(np.abs(dl)>1e-6)
            if sub.any():
                ss=np.abs(dv[sub]/dl[sub])*1e3; SS=float(ss.min()); par.append(f"Subthreshold swing SS = {SS:.1f} mV/dec")
            gm=np.gradient(Ia,Vg); j=int(np.argmax(np.abs(gm)))
            if gm[j]!=0:
                Vth=Vg[j]-Ia[j]/gm[j]
                edge=j in (0,len(Vg)-1)
                par.append(f"V_th (max-g_m linear extrapolation) = {Vth:.3f} V"+
                           ("  - g_m still rising at the sweep end: extend the sweep for a firm value" if edge else ""))
                par.append(f"g_m,max = {abs(gm[j]):.3e} S per cell at V = {Vg[j]:.3f} V")
                aL.axvline(Vth,color=C["yellow"],lw=1,ls="--",alpha=0.7)
                aL.text(Vth,aL.get_ylim()[1],f" V_th {Vth:.2f} V",color=C["yellow"],fontsize=8,va="top")
            if wl:
                Icc=1e-7*wl; k=np.where((Ia[:-1]<Icc)&(Ia[1:]>=Icc)|(Ia[:-1]>=Icc)&(Ia[1:]<Icc))[0]
                if len(k):
                    k=k[0]; l0,l1=np.log(Ia[k]),np.log(Ia[k+1])
                    Vcc=Vg[k]+(Vg[k+1]-Vg[k])*(np.log(Icc)-l0)/(l1-l0)
                    par.append(f"V_th (constant current 100 nA x W/L = {Icc:.2e} A) = {Vcc:.3f} V")
            Ion,Ioff=Ia[-1],Ia[0]
            if abs(Vg[0])>abs(Vg[-1]): Ion,Ioff=Ioff,Ion
            par.append(f"I_on = {Ion:.3e} A, I_off = {Ioff:.3e} A per cell"+(f"  (I_on/I_off = {Ion/Ioff:.2e})" if Ioff>0 else ""))
            if par: aS.text(0.03,0.97,"\n".join(p.split("  - ")[0] for p in par[:2]),transform=aS.transAxes,va="top",fontsize=8,color=C["yellow"],
                            bbox=dict(boxstyle="round,pad=0.4",facecolor=C["panel"],edgecolor=C["border"],alpha=0.9))
        Aok_cm2=Aok*1e4
        if tkey==okey and Aok_cm2>0 and not gate.get(okey,False):
            J=Io/Aok_cm2
            fwd=good&(J>1e-20)&(V>0.05)
            if fwd.sum()>=3 and np.log10(J[fwd].max()/J[fwd].min())>=2.:
                try:
                    m_fit,b=np.polyfit(V[fwd],np.log(J[fwd]),1); n_id=1./(m_fit*0.02585); J0=np.exp(b)
                    if 0.8<n_id<4.:
                        par.append(f"Diode ideality n = {n_id:.3f}, J0 = {J0:.2e} A/cm²")
                        aS.text(0.03,0.97,f"n = {n_id:.3f}\nJ₀ = {J0:.2e} A/cm²",transform=aS.transAxes,va="top",
                                fontsize=8,color=C["yellow"],bbox=dict(boxstyle="round,pad=0.4",facecolor=C["panel"],
                                edgecolor=C["border"],alpha=0.9))
                except Exception: pass
        if bv>0 and Aok_cm2>0 and np.max(np.abs(V))>=3.:
            Ja=np.where(good,np.abs(Io)/Aok_cm2,0.); Vbd=None
            for k in range(1,len(Ja)):
                if Ja[k-1]>1e-30 and Ja[k]/Ja[k-1]>10. and Ja[k]>bv: Vbd=V[k]; break
            if Vbd is None:
                for k in range(len(Ja)):
                    if Ja[k]>bv: Vbd=V[k]; break
            if Vbd is not None:
                for ax in (aL,aS): ax.axvline(Vbd,color="#ff4444",lw=2,ls="--",alpha=0.8)
                par.append(f"Current knee at V = {Vbd:.2f} V (|J| > {bv:g} A/cm²; heuristic - no avalanche model)")

    def _plot_vtc(self,fig,V,Vf,cv,gI,gF,gate,tkey,fkey,okey,par):
        gs=GridSpec(1,2,figure=fig,wspace=0.45)
        aT=fig.add_subplot(gs[0]); aI=fig.add_subplot(gs[1])
        aT.set_facecolor(C["bg"]); aI.set_facecolor(C["bg"])
        ok=cv&np.isfinite(Vf)
        Vi,idx=np.unique(V[ok],return_index=True); Vo=Vf[ok][idx]
        aT.plot(V[ok],Vf[ok],color=C["green"],lw=2,marker="o",ms=3,label=f"V({fkey})")
        if (~cv).any(): aT.plot(V[~cv],Vf[~cv],color=C["red"],marker="x",ls="none",ms=6,label="not converged")
        if len(V): aT.plot([V.min(),V.max()],[V.min(),V.max()],color=C["sub"],ls=":",lw=1,label="V_out = V_in")
        if len(Vi)>=3:
            g=np.gradient(Vo,Vi); a2=aT.twinx()
            a2.plot(Vi,g,color=C["orange"],lw=1.2,ls="--",label="gain")
            a2.set_ylabel("gain dV_out/dV_in",fontsize=8,color=C["orange"]); a2.tick_params(labelsize=7,colors=C["orange"])
            d_=Vo-Vi; k=np.where(np.diff(np.sign(d_))!=0)[0]
            if len(k):
                k=k[0]; VM=Vi[k]+(Vi[k+1]-Vi[k])*d_[k]/(d_[k]-d_[k+1])
                aT.axvline(VM,color=C["yellow"],lw=1,ls="--",alpha=0.8)
                aT.text(VM,max(Vo),f" V_M = {VM:.3f} V",color=C["yellow"],fontsize=8,va="top")
                par.append(f"Switching threshold V_M = {VM:.3f} V")
            gi=int(np.argmax(np.abs(g))); par.append(f"Max |gain| = {abs(g[gi]):.1f} at V_in = {Vi[gi]:.3f} V")
            VOH,VOL=(Vo[0],Vo[-1]) if Vo[0]>Vo[-1] else (Vo[-1],Vo[0])
            par.append(f"V_OH = {VOH:.4f} V, V_OL = {VOL:.4f} V")
            c=np.where(np.diff(np.sign(g+1.))!=0)[0]
            if len(c)>=2:
                f=lambda j:Vi[j]+(Vi[j+1]-Vi[j])*(-1.-g[j])/(g[j+1]-g[j])
                VIL,VIH=f(c[0]),f(c[-1])
                par.append(f"V_IL = {VIL:.3f} V, V_IH = {VIH:.3f} V (gain = -1)")
                par.append(f"Noise margins NM_L = {VIL-VOL:.3f} V, NM_H = {VOH-VIH:.3f} V")
                for vv in (VIL,VIH): aT.axvline(vv,color=C["sub"],lw=0.8,ls=":")
        aT.legend(fontsize=7,loc="center left")
        sax(aT,f"V({tkey}) (V)",f"V({fkey}) (V)","Voltage transfer curve (output floating, I = 0)")
        for n,k in enumerate([k for k in gI if not gate[k] and k!=fkey]):
            Ic=np.abs(gI[k]); okk=cv&(Ic>=gF[k])&(Ic>0)
            if okk.any():
                aI.semilogy(V[okk],Ic[okk],color=TCOLS[n%len(TCOLS)],lw=2.2 if k==okey else 1.2,
                            marker="o" if k==okey else None,ms=3,label=f"|I| {k}")
        if np.any(gF[okey]>0): aI.semilogy(V,np.maximum(gF[okey],1e-300),color=C["sub"],lw=0.8,ls=":",label="resolution")
        aI.legend(fontsize=7); sax(aI,f"V({tkey}) (V)","|I| (A per cell)","Supply / terminal currents")
        Ia=np.abs(gI[okey])[cv]
        if len(Ia):
            j=int(np.argmax(Ia)); par.append(f"Peak |I({okey})| = {Ia[j]:.3e} A per cell at V_in = {V[cv][j]:.3f} V")

    # ─── EXPORT / FILES ────────────────────────────────────────
    def _xcsv(self):
        s=self.sol
        if s is None: messagebox.showinfo("Export","Run a simulation first."); return
        p=filedialog.asksaveasfilename(defaultextension=".csv",filetypes=[("CSV","*.csv")],
                                       initialfile=f"{self.tname or 'device'}_solution.csv".replace(" ","_").replace(":",""))
        if not p: return
        keys=("phi","n","p","Ec","Ev","Efn","Efp","Jx","Jy","Jz","mn","mp","R","Nnet")
        if s.get("mesh")=="tri" and "raw" in s:        # the prism mesh's own nodes
            P=s["prism"].node_xyz(); cols=[P[:,0],P[:,1],P[:,2]]+[s["raw"][k] for k in keys]
        else:
            Z,Y,X=np.meshgrid(s["z"],s["y"],s["x"],indexing="ij")
            cols=[X,Y,Z]+[s[k] for k in keys]
        data=np.column_stack([np.asarray(c).reshape(-1) for c in cols])
        np.savetxt(p,data,delimiter=",",fmt="%.6e",comments="",
                   header="x_nm,y_nm,z_nm,phi_V,n_cm3,p_cm3,Ec_eV,Ev_eV,Efn_eV,Efp_eV,Jx_Acm2,Jy_Acm2,Jz_Acm2,"
                          "mun_cm2Vs,mup_cm2Vs,R_cm3s,Nnet_cm3")
        self.sv.set(f"Solution → {p}"); self._log(f"Solution exported to {p}")
    def _xiv_csv(self):
        d=self.ivd
        if d is None: messagebox.showinfo("Export","Run a sweep first."); return
        p=filedialog.asksaveasfilename(defaultextension=".csv",filetypes=[("CSV","*.csv")],
                                       initialfile=f"{self.tname or 'device'}_sweep.csv".replace(" ","_").replace(":",""))
        if not p: return
        R=d["R"]; labs=d["labels"]
        cols=[R["V"]]; hdr=[f"V_{d['tkey']}_V"]
        if d["fkey"]: cols.append(R["Vf"]); hdr.append(f"V_{d['fkey']}_V")
        cols.append(R["conv"].astype(float)); hdr.append("converged")
        for c,l in enumerate(labs):
            cols.append(R["I"][:,c]); hdr.append(f"I_{l}_A")
        for c,l in enumerate(labs):
            cols.append(R["F"][:,c]); hdr.append(f"floor_{l}_A")
        np.savetxt(p,np.column_stack(cols),delimiter=",",fmt="%.8e",comments="",
                   header=",".join(h.replace(",",";") for h in hdr))
        self.sv.set(f"Sweep → {p}"); self._log(f"Sweep exported to {p}")
    def _xplt(self):
        p=filedialog.asksaveasfilename(defaultextension=".png",filetypes=[("PNG","*.png"),("PDF","*.pdf"),("SVG","*.svg")])
        if not p: return
        figs={str(self.td):self.fd,str(self.ts):self.fsl,str(self.tb):self.fb,str(self.tc):self.fc,
              str(self.ti):self.fiv,str(self.tf):self.ff}
        fig=figs.get(self.nb.select())
        if fig is None: messagebox.showinfo("Export","The visible tab has no plot."); return
        fig.savefig(p,dpi=200,bbox_inches="tight",facecolor=C["bg"])
        self.sv.set(f"Plot → {p}")

    def _save_device(self):
        p=filedialog.asksaveasfilename(defaultextension=".json",filetypes=[("SemiSim device","*.json")],
                                       initialfile=f"{self.tname or 'device'}.json".replace(" ","_").replace(":",""))
        if not p: return
        try:
            d=dict(format="semisim3d-device",version=1,name=self.tname,
                   Lx=self.Lx.get(),Ly=self.Ly.get(),Lz=self.Lz.get(),Nx=self.Nx.get(),Ny=self.Ny.get(),Nz=self.Nz.get(),
                   graded=self.use_graded.get(),dense_x=self.dense_x,dense_y=self.dense_y,dense_z=self.dense_z,
                   meta=self.tmeta,regs=[obj_dict(r) for r in self.regs],cons=[obj_dict(c) for c in self.cons],
                   sheets=[obj_dict(s) for s in self.sheets],
                   mesh=dict(type=self.mesh_type.get(),prism_axis=self.prism_ax.get()),
                   sweep=dict(tgt=self._gmap.get(self.iv_tgt.get()),flt=self._gmap.get(self.iv_flt.get()),
                              out=self._gmap.get(self.iv_out.get()),Vs=self.ivS.get(),Ve=self.ivE.get(),
                              N=self.ivN.get(),bv=self.bv_thr.get()),
                   solver=dict(T=self.T.get(),max_p=self.max_p.get(),max_c=self.max_c.get(),max_g=self.max_g.get(),
                               tol_p=self.tol_p.get(),tol_c=self.tol_c.get(),bgn=self.m_bgn.get(),
                               tau=self.m_tau.get(),vsat=self.m_vsat.get(),qc=self.m_qc.get(),
                               lombardi=self.m_lb.get(),tunnel=self.m_tun.get()))
            with open(p,"w") as f: json.dump(d,f,indent=1)
            self.sv.set(f"Device saved → {p}"); self._log(f"Device saved to {p}")
        except Exception:
            messagebox.showerror("Save device",traceback.format_exc())
    def _open_device(self):
        if self.busy: return
        p=filedialog.askopenfilename(filetypes=[("SemiSim device","*.json"),("All files","*.*")])
        if not p: return
        try:
            with open(p) as f: d=json.load(f)
            if d.get("format")!="semisim3d-device": raise ValueError("not a SemiSim 3D device file")
            regs=[Reg(**r) for r in d["regs"]]; cons=[Con(**c) for c in d["cons"]]
            sheets=[Sheet(**s) for s in d.get("sheets",[])]
        except Exception as e:
            messagebox.showerror("Open device",f"Cannot read {p}:\n{e}"); return
        self.tname=d.get("name") or os.path.basename(p); self.tmeta=d.get("meta") or {}
        for v,k in ((self.Lx,"Lx"),(self.Ly,"Ly"),(self.Lz,"Lz"),(self.Nx,"Nx"),(self.Ny,"Ny"),(self.Nz,"Nz")): v.set(d[k])
        self.use_graded.set(d.get("graded",True))
        ms=d.get("mesh") or {}
        self.mesh_type.set(ms.get("type","rect")); self.prism_ax.set(ms.get("prism_axis","auto"))
        self.dense_x,self.dense_y,self.dense_z=d.get("dense_x",[]),d.get("dense_y",[]),d.get("dense_z",[])
        self.regs,self.cons,self.sheets=regs,cons,sheets
        sw=d.get("sweep",{})
        self.ivS.set(sw.get("Vs",0.)); self.ivE.set(sw.get("Ve",1.)); self.ivN.set(sw.get("N",11)); self.bv_thr.set(sw.get("bv",0.))
        so=d.get("solver",{})
        for v,k in ((self.T,"T"),(self.max_p,"max_p"),(self.max_c,"max_c"),(self.max_g,"max_g"),
                    (self.tol_p,"tol_p"),(self.tol_c,"tol_c"),(self.m_bgn,"bgn"),(self.m_tau,"tau"),(self.m_vsat,"vsat"),
                    (self.m_qc,"qc"),(self.m_lb,"lombardi"),(self.m_tun,"tunnel")):
            v.set(so[k] if k in so else (False if k in ("qc","lombardi","tunnel") else v.get()))
        self._refresh_groups(tgt=sw.get("tgt"),flt=sw.get("flt") or "(none)",out=sw.get("out") or "auto")
        self._rfr(redraw=False); self._rfc(redraw=False); self._show_sheets()
        self.t_name.config(text=self.tname); self.t_cat.config(text=f"from {os.path.basename(p)}")
        self.t_doc.configure(state="normal"); self.t_doc.delete("1.0","end")
        self.t_doc.insert("end",template_doc(self.tname) if self.tname in TEMPLATES else f"Device loaded from {p}")
        self.t_doc.configure(state="disabled"); self.tmpl_btn.config(text=self.tname[:30])
        self.sol=None; self.ivd=None; self._clear_solution_plots(); self._fill_results(None)
        self._fresh_cuts=True
        view=self.tmeta.get("view")
        if view: self.sec_plane.set(view[0]); self._sec_range(view[1])
        else: self._sec_range(None)
        self._mesh_update(); self.sv.set(f"Opened {p}"); self._log(f"Device opened from {p}")

    # ─── HELP / ABOUT / TEMPLATE BROWSER ───────────────────────
    def _text_window(self,title,w=900,h=640):
        top=tk.Toplevel(self); top.title(title); top.configure(bg=C["panel"])
        sw,sh=self.winfo_screenwidth(),self.winfo_screenheight()
        W,H=min(w,sw-60),min(h,sh-100)
        top.geometry(f"{W}x{H}+{self.winfo_rootx()+40}+{self.winfo_rooty()+30}")
        top.minsize(420,300); top.transient(self)
        top.bind("<Escape>",lambda e:top.destroy())
        return top
    def _styled_text(self,parent):
        fr=tk.Frame(parent,bg=C["panel"])
        t=tk.Text(fr,wrap="word",bg="#07111e",fg=C["text"],relief="flat",font=self.F["ui"],padx=14,pady=10,
                  highlightthickness=0,spacing1=2,spacing3=2,insertbackground=C["text"])
        vs=ttk.Scrollbar(fr,orient="vertical",command=t.yview); t.configure(yscrollcommand=vs.set)
        t.pack(side=tk.LEFT,fill=tk.BOTH,expand=True); vs.pack(side=tk.RIGHT,fill=tk.Y)
        t.tag_configure("h1",font=self.F["big"],foreground=C["accent"],spacing1=10,spacing3=6)
        t.tag_configure("h2",font=self.F["head"],foreground=C["yellow"],spacing1=8,spacing3=3)
        t.tag_configure("bullet",lmargin1=10,lmargin2=30)
        t.tag_configure("bullet2",lmargin1=10,lmargin2=52)
        t.tag_configure("code",font=self.F["mono"],foreground="#86efac")
        t.tag_configure("pre",font=self.F["mono"],foreground="#cbd5e1",lmargin1=16,lmargin2=16)
        t.tag_configure("hit",background="#6b5a00",foreground="#ffffff")
        return fr,t

    def _help(self,topic=None):
        if self._helpwin is not None and self._helpwin.winfo_exists():
            self._helpwin.deiconify(); self._helpwin.lift()
            if topic: self._help_goto(topic)
            return
        top=self._text_window(f"SemiSim 3D {VERSION} - Help",980,700); self._helpwin=top
        bar=tk.Frame(top,bg=C["panel"]); bar.pack(fill=tk.X,padx=6,pady=4)
        tk.Label(bar,text="Find",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=4)
        q=tk.StringVar(); e=tk.Entry(bar,textvariable=q,bg="#07111e",fg=C["text"],insertbackground=C["text"],
                                     relief="flat",font=self.F["ui"],width=28)
        e.pack(side=tk.LEFT,padx=4)
        pw=ttk.PanedWindow(top,orient=tk.HORIZONTAL); pw.pack(fill=tk.BOTH,expand=True,padx=6,pady=(0,6))
        lf=tk.Frame(pw,bg=C["panel"])
        lb=tk.Listbox(lf,bg="#07111e",fg=C["text"],selectbackground=C["accent"],relief="flat",font=self.F["ui"],
                      activestyle="none",highlightthickness=0,width=24)
        lb.pack(fill=tk.BOTH,expand=True)
        fr,t=self._styled_text(pw)
        pw.add(lf,weight=0); pw.add(fr,weight=1)
        self._help_marks=[]
        for k,(title,body) in enumerate(help_topics()):
            mk=f"topic{k}"; t.mark_set(mk,"end-1c"); t.mark_gravity(mk,"left")
            self._help_marks.append((title,mk)); lb.insert("end",title)
            render_markup(t,body); t.insert("end","\n")
        t.configure(state="disabled"); self._help_text=t; self._help_list=lb
        lb.bind("<<ListboxSelect>>",lambda _:(lb.curselection() and t.yview(self._help_marks[lb.curselection()[0]][1])))
        state={"last":"1.0"}
        def find(_=None):
            t.tag_remove("hit","1.0","end"); s=q.get().strip()
            if not s: return
            idx="1.0"; first=None
            while True:
                idx=t.search(s,idx,nocase=True,stopindex="end")
                if not idx: break
                end=f"{idx}+{len(s)}c"; t.tag_add("hit",idx,end); first=first or idx; idx=end
            nxt=t.search(s,state["last"],nocase=True,stopindex="end") or first
            if nxt: t.see(nxt); state["last"]=f"{nxt}+{len(s)}c"
        e.bind("<Return>",find)
        self._btn(bar,"Find next",find,bg=C["border"]).pack(side=tk.LEFT,padx=4)
        self._btn(bar,"Close",top.destroy,bg=C["border"]).pack(side=tk.RIGHT,padx=4)
        if topic: self._help_goto(topic)
    def _help_goto(self,topic):
        for k,(title,mk) in enumerate(self._help_marks):
            if title==topic:
                self._help_text.yview(mk); self._help_list.selection_clear(0,"end"); self._help_list.selection_set(k)

    def _about(self):
        top=self._text_window(f"About SemiSim 3D {VERSION}",720,560)
        fr,t=self._styled_text(top); fr.pack(fill=tk.BOTH,expand=True,padx=6,pady=6)
        body=f"""
# SemiSim 3D {VERSION}
3-D drift-diffusion semiconductor device simulator: a C core (OpenMP) with a Python/Tk interface.
## This session
- core library: `{_SO}`
- OpenMP threads: {os.environ.get('OMP_NUM_THREADS','?')} (of {self.ncpu} available)
- memory: about {BYTES_PER_NODE} bytes per mesh node; Anderson depth {api_aa() or '-'}
- settings file: `{SETTINGS_FILE}`
## Transport core (v7)
- box method on rectangular (tensor) meshes or triangular-prism meshes (Delaunay cross-section x 1-D grid)
- Newton-Poisson (multigrid CG) + Scharfetter-Gummel continuity (multigrid BiCGSTAB)
- Gummel map on (ψ, φn, φp) with Anderson mixing, floating-region balance
- adaptive bias continuation, cooperative stop, Kirchhoff-exact terminal currents with resolution floor
## Physics
- heterojunctions, oxide gates (tox, Nox, layer stacks), SiO2/Al2O3/HfO2/Si3N4, buried electrodes, interface sheet charge
- Arora mobility + velocity saturation + Lombardi surface mobility, SRH τ(N) + Auger + radiative, Klaassen band-gap narrowing
- quantum confinement (MLDA), gate tunnelling (Tsu-Esaki/WKB: direct + Fowler-Nordheim)
## Not modelled
- impact ionisation, band-to-band tunnelling, optical generation, subband quantisation
- thermionic heterointerfaces, ballistic transport, self-heating
## Build
  gcc -O3 -march=native -ffast-math -fno-finite-math-only -fopenmp -shared -fPIC \\
      -o semiconductor_core.so semiconductor_3d.c -lm
- `semimesh.py` (the triangular-prism mesher) must sit next to main.py
"""
        render_markup(t,body); t.configure(state="disabled")

    def _browse_templates(self):
        if self._browser is not None and self._browser.winfo_exists():
            self._browser.deiconify(); self._browser.lift(); return
        top=self._text_window("Device template library",1000,640); self._browser=top
        pw=ttk.PanedWindow(top,orient=tk.HORIZONTAL); pw.pack(fill=tk.BOTH,expand=True,padx=6,pady=6)
        lf=tk.Frame(pw,bg=C["panel"])
        tv=ttk.Treeview(lf,show="tree",selectmode="browse")
        vs=ttk.Scrollbar(lf,orient="vertical",command=tv.yview); tv.configure(yscrollcommand=vs.set)
        tv.pack(side=tk.LEFT,fill=tk.BOTH,expand=True); vs.pack(side=tk.RIGHT,fill=tk.Y)
        tv.column("#0",width=300)
        for cat,names in template_groups():
            p=tv.insert("","end",text=cat,open=True)
            for n in names: tv.insert(p,"end",text=n,values=(n,))
        rf=tk.Frame(pw,bg=C["panel"])
        fr,t=self._styled_text(rf); fr.pack(fill=tk.BOTH,expand=True)
        bb=tk.Frame(rf,bg=C["panel"]); bb.pack(fill=tk.X,pady=4)
        sel={"n":None}
        def show(_=None):
            it=tv.selection()
            if not it: return
            v=tv.item(it[0],"values")
            if not v: return
            n=v[0]; sel["n"]=n
            t.configure(state="normal"); t.delete("1.0","end")
            t.insert("end",n+"\n","h1"); t.insert("end",template_doc(n)); t.configure(state="disabled")
        def load(_=None):
            if sel["n"]: self._load(sel["n"])
        tv.bind("<<TreeviewSelect>>",show); tv.bind("<Double-1>",load)
        self._btn(bb,"Load this device",load,bg="#0a8f5a",font=self.F["uib"]).pack(side=tk.LEFT,padx=6)
        self._btn(bb,"Close",top.destroy,bg=C["border"]).pack(side=tk.RIGHT,padx=6)
        pw.add(lf,weight=0); pw.add(rf,weight=1)


if __name__=="__main__":
    if sys.platform.startswith("win"):
        try: ctypes.windll.shcore.SetProcessDpiAwareness(1)      # crisp text on high-DPI screens
        except Exception: pass
    App().mainloop()
