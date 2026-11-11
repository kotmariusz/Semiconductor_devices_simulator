#!/usr/bin/env python3
"""SemiSim 3D — 3-D drift-diffusion semiconductor device simulator.
Rectangular (graded tensor) or triangular-prism meshes (semimesh.py);
dense nodes at junctions, interfaces and gates, coarse elsewhere."""

import sys,os,ctypes,threading,time,traceback,platform,subprocess,math
IS_MAC=sys.platform=="darwin"
def _is_wsl():
    """Windows Subsystem for Linux (its WSLg window system has quirks)."""
    if not sys.platform.startswith("linux"): return False
    if os.environ.get("WSL_DISTRO_NAME") or os.environ.get("WSL_INTEROP"): return True
    try: return "microsoft" in platform.release().lower()
    except Exception: return False
IS_WSL=_is_wsl()
def _mac_perf_cores():
    """Performance cores of an Apple-silicon Mac (0 if unknown).  The solver's
    parallel loops wait for the slowest thread at every barrier, so the
    efficiency cores slow it down rather than help."""
    try:
        r=subprocess.run(["sysctl","-n","hw.perflevel0.physicalcpu"],capture_output=True,text=True,timeout=5)
        return int(r.stdout.strip() or 0)
    except Exception: return 0
# Use all available CPU cores for OpenMP-parallel C solver, unless user
# already set OMP_NUM_THREADS in their environment.
if "OMP_NUM_THREADS" not in os.environ:
    try:
        n_cpu = len(os.sched_getaffinity(0))   # respects taskset/cpusets
    except (AttributeError, OSError):
        n_cpu = os.cpu_count() or 1
    if IS_MAC: n_cpu = _mac_perf_cores() or n_cpu
    os.environ["OMP_NUM_THREADS"] = str(max(1, n_cpu))
# Short OpenMP spin-wait before sleeping: ~3% slower on an idle machine, but
# 2-3x faster when other programs compete for the cores (long spinning
# threads then wait on descheduled ones at every barrier).
os.environ.setdefault("GOMP_SPINCOUNT", "30000")

_APP_NAP_TOKEN=None
def _mac_no_app_nap():
    """macOS: opt out of App Nap, which would throttle the solver (and move
    it to the efficiency cores) while the window is hidden behind others.
    [[NSProcessInfo processInfo] beginActivityWithOptions:reason:] through
    the Objective-C runtime, as the 'appnope' package does."""
    global _APP_NAP_TOKEN
    if not IS_MAC or _APP_NAP_TOKEN is not None: return
    try:
        import ctypes.util
        objc=ctypes.cdll.LoadLibrary(ctypes.util.find_library("objc"))
        cf=ctypes.cdll.LoadLibrary(ctypes.util.find_library("CoreFoundation"))
        ctypes.cdll.LoadLibrary(ctypes.util.find_library("Foundation"))
        vp=ctypes.c_void_p
        objc.objc_getClass.restype=vp; objc.objc_getClass.argtypes=[ctypes.c_char_p]
        objc.sel_registerName.restype=vp; objc.sel_registerName.argtypes=[ctypes.c_char_p]
        msg0=ctypes.CFUNCTYPE(vp,vp,vp)(("objc_msgSend",objc))
        msg2=ctypes.CFUNCTYPE(vp,vp,vp,ctypes.c_uint64,vp)(("objc_msgSend",objc))
        cfstr=cf.CFStringCreateWithCString; cfstr.restype=vp; cfstr.argtypes=[vp,ctypes.c_char_p,ctypes.c_uint32]
        pinfo=msg0(objc.objc_getClass(b"NSProcessInfo"),objc.sel_registerName(b"processInfo"))
        opts=0x00FFFFFF & ~(1<<20)            # NSActivityUserInitiatedAllowingIdleSystemSleep
        _APP_NAP_TOKEN=msg2(pinfo,objc.sel_registerName(b"beginActivityWithOptions:reason:"),opts,
                            cfstr(None,b"SemiSim 3D device simulation",0x08000100))
    except Exception:
        _APP_NAP_TOKEN=0

def _fatal(msg):
    """Report a start-up problem in the terminal and, when possible, in a window."""
    print("\nSemiSim 3D cannot start:\n"+msg+"\n",file=sys.stderr)
    try:
        import tkinter as _tk; from tkinter import messagebox as _mb
        r=_tk.Tk(); r.withdraw(); _mb.showerror("SemiSim 3D cannot start",msg); r.destroy()
    except Exception: pass
    sys.exit(1)

try:
    import tkinter as tk
    from tkinter import ttk,messagebox,filedialog
except ImportError:
    _fatal("Python's Tk interface (tkinter) is missing.\n"
           +("Install Python from python.org (it includes Tk), or:  brew install python-tk"
             if IS_MAC else "Ubuntu / Debian / WSL:  sudo apt install python3-tk"))
try:
    import numpy as np
    import matplotlib
except ImportError as e:
    _fatal(f"A Python package is missing ({e}).\nInstall everything with:  bash install.sh\n"
           "or:  python3 -m pip install -r requirements.txt")
if tk.TkVersion<8.6:
    _fatal(f"This Python uses Tcl/Tk {tk.TkVersion}, which is too old (macOS's built-in Python).\n"
           "Install Python 3 from python.org (or Homebrew) and run install.sh again.")
def _ver(s):
    out=[]
    for p in s.split(".")[:2]:
        d="".join(ch for ch in p if ch.isdigit()); out.append(int(d or 0))
    return tuple(out)
if tk.TkVersion>=9.0 and _ver(matplotlib.__version__)<(3,10):
    _fatal(f"Tcl/Tk {tk.TkVersion} needs matplotlib 3.10 or newer (found {matplotlib.__version__}).\n"
           "Update it:  python3 -m pip install -U matplotlib")
matplotlib.use("TkAgg")
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
    import semidefects as SD               # defects: library, shapes, mapping onto a mesh
except ImportError as e:
    _fatal(f"semimesh.py and semidefects.py must be in the same folder as main.py ({e})")

# ══════════════════════════════════════════════════════════════
#  C LIBRARY
# ══════════════════════════════════════════════════════════════
import semibuild                          # builds the core on the first start (same folder)
try: _lib,_SO=semibuild.load_core()
except semibuild.CoreError as e: _fatal(str(e))

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
api_conIraw = _f("api3_get_contact_Iraw",D,I)   # the contact's own flux (A), before the Kirchhoff choice
api_stat    = _f("api3_get_stat",ctypes.c_long,I)
api_Nx=_f("api3_get_Nx",I); api_Ny=_f("api3_get_Ny",I); api_Nz=_f("api3_get_Nz",I)
api_conv=_f("api3_get_conv",I); api_iter=_f("api3_get_iter",I); api_res=_f("api3_get_res",D)
api_xs=_f("api3_get_xs",D,I); api_ys=_f("api3_get_ys",D,I); api_zs=_f("api3_get_zs",D,I)
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
BYTES_PER_NODE=760   # the core incl. multigrid hierarchy and depth-12 Anderson history (depth drops if RAM is short)

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
         "SiO2":8,"Al2O3":9,"HfO2":10,"Si3N4":11,"Void":12,"Metal (incl.)":13}
INSULATORS = {"SiO2","Al2O3","HfO2","Si3N4","Void","Metal (incl.)"}
EPS_R = {"SiO2":3.9,"Al2O3":9.0,"HfO2":22.0,"Si3N4":7.5,"Void":1.0,"Metal (incl.)":1000.}   # insulators: EOT of gate stacks

# ── contact metals ─────────────────────────────────────────────
# Work function (eV, polycrystalline films; Michaelson 1977) and the measured
# Schottky barrier heights on n-type semiconductors that are well established
# (eV; Si, GaAs and silicides: Sze & Ng, Physics of Semiconductor Devices;
# 4H-SiC and GaN: typical values of the contact literature).  Any other pair
# uses the Cowley-Sze pinning model below.  Real barriers depend on surface
# preparation and anneal (+-0.1 eV): pick "(custom phi_m)" and type phi_m to
# override.  Gates use the work function (poly-Si: its band edge).
METALS = {
    "Al":   (4.28, {"Si":0.72, "GaAs":0.80, "4H-SiC":1.00, "GaN":0.60}),
    "Ti":   (4.33, {"Si":0.50, "GaAs":0.83, "4H-SiC":1.10, "GaN":0.58}),
    "Ag":   (4.26, {"Si":0.78, "GaAs":0.88}),
    "Cr":   (4.50, {"Si":0.61}),
    "W":    (4.55, {"Si":0.67, "GaAs":0.80, "4H-SiC":1.15}),
    "Mo":   (4.60, {"Si":0.68, "4H-SiC":1.20}),
    "Cu":   (4.65, {"Si":0.58, "GaAs":0.82}),
    "Au":   (5.10, {"Si":0.80, "GaAs":0.90, "4H-SiC":1.60, "GaN":0.95}),
    "Pd":   (5.12, {"Si":0.81, "GaN":0.95}),
    "Ni":   (5.15, {"Si":0.61, "GaAs":0.83, "4H-SiC":1.50, "GaN":1.00}),
    "Pt":   (5.65, {"Si":0.90, "GaAs":0.86, "4H-SiC":1.60, "GaN":1.10}),
    "TiN":  (4.60, {}),
    "TaN":  (4.50, {}),
    "TiSi2":(4.60, {"Si":0.60}),
    "CoSi2":(4.70, {"Si":0.64}),
    "NiSi": (4.80, {"Si":0.65}),
    "PtSi": (4.95, {"Si":0.84}),
    "ErSi2":(3.90, {"Si":0.28}),
    "n+ poly-Si":(4.05, {}),
    "p+ poly-Si":(5.17, {}),
}
METAL_NONE="(ideal / custom φm)"
# semiconductor electron affinity and gap (eV, 300 K) as in the core
SEMI_CHI = {"Si":4.05,"Ge":4.00,"GaAs":4.07,"InP":4.38,"4H-SiC":3.70,"GaN":4.10,"AlGaAs":3.827,"AlGaN":3.76}
SEMI_EG  = {"Si":1.124,"Ge":0.66,"GaAs":1.424,"InP":1.353,"4H-SiC":3.23,"GaN":3.43,"AlGaAs":1.80,"AlGaN":3.90}
# Cowley-Sze: phi_Bn = S (phi_m - chi) + (1 - S) phi_pin  (slope S, pinned barrier phi_pin)
SB_PIN   = {"Si":(0.27,0.80),"Ge":(0.05,0.58),"GaAs":(0.07,0.85),"InP":(0.10,0.48),
            "4H-SiC":(0.40,1.25),"GaN":(0.40,0.90),"AlGaAs":(0.07,1.00),"AlGaN":(0.40,1.20)}

def schottky_barrier(metal, semi):
    """(phi_Bn in eV, source) of a metal on the n-type semiconductor `semi`;
    the hole barrier is Eg - phi_Bn."""
    wf,meas=METALS[metal]
    if semi in meas: return meas[semi],"measured"
    S,pin=SB_PIN.get(semi,(1.,0.))
    phi=S*(wf-SEMI_CHI.get(semi,4.05))+(1.-S)*pin
    return float(min(max(phi,0.),SEMI_EG.get(semi,1.12))),"pinning model"

def contact_semiconductor(con, regs, L):
    """Material of the region under a contact (the middle of its face, just
    inside the device; a box: its centre); None if an insulator or no region."""
    if not regs: return None
    e=con_extent(con,*L); P={a:0.5*(e[a][0]+e[a][1]) for a in "xyz"}
    if con.face<=5:
        ai=con.face//2; P["xyz"[ai]] = 1e-4*L[ai] if con.face%2==0 else L[ai]*(1.-1e-4)
    idx=int(np.asarray(SM.region_index(regs,np.array([P["x"]]),np.array([P["y"]]),np.array([P["z"]]),L)).ravel()[0])
    if idx<0: return None
    m=regs[idx].mat
    return None if m in INSULATORS else m

def contact_core(con, regs, L, models=None):
    """(bc, phi_m) the core gets for a contact, and a description.  A contact
    with a metal: a gate takes the metal's work function; on a semiconductor
    (Ohmic or Schottky) it becomes a physical metal contact (core bc 1) whose
    phi_m puts the Schottky barrier of that metal at the interface - thermionic
    emission and barrier tunnelling then decide whether it behaves ohmic.
    With the metal-contact model off (models without bit 64) an Ohmic contact
    stays ideal and a Schottky contact is the classic Dirichlet barrier."""
    metal=getattr(con,"metal","") or ""
    if metal not in METALS: return con.bc, con.pm, ""
    wf=METALS[metal][0]
    if con.bc==2: return 2, wf, f"{metal}: φm {wf:.2f} eV"
    semi=contact_semiconductor(con,regs,L)
    if semi is None or semi not in SEMI_CHI: return con.bc, wf, f"{metal}: φm {wf:.2f} eV (not on a semiconductor)"
    pb,src=schottky_barrier(metal,semi)
    txt=f"{metal} on {semi}: φBn {pb:.2f} eV, φBp {SEMI_EG[semi]-pb:.2f} eV ({src})"
    if con.face>5:
        # buried (box) contact: no emission/tunnelling model - an Ohmic box
        # stays ideal, a Schottky box is the fixed-density barrier of the metal
        if con.bc==0: return 0, con.pm, txt+" - box contact: ideal ohmic"
        return 1, SEMI_CHI[semi]+pb, txt+" - box contact: fixed-density barrier"
    if con.bc==0 and models is not None and not (int(models)&64):
        return 0, con.pm, txt+" - ideal ohmic while the metal-contact model is off"
    return 1, SEMI_CHI[semi]+pb, txt
MCOLS = {"Si":("#3b9de8","#e85050"),"Ge":("#2ec08a","#f07030"),
         "GaAs":("#a855f7","#e8305a"),"InP":("#06b6d4","#e05090"),
         "4H-SiC":("#22c55e","#eab308"),"GaN":("#38bdf8","#f97316"),
         "AlGaAs":("#8b5cf6","#f43f5e"),"AlGaN":("#67e8f9","#fb923c"),
         "SiO2":("#cbd5e1","#cbd5e1"),"Al2O3":("#e2e8f0","#e2e8f0"),
         "HfO2":("#fcd34d","#fcd34d"),"Si3N4":("#d6d3d1","#d6d3d1"),
         "Void":("#111827","#111827"),"Metal (incl.)":("#d4a017","#d4a017")}
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

def _cur_unit(m):
    """(scale, unit) that puts a current of magnitude m (A) between 1 and 1000."""
    for sc,un in ((1.,"A"),(1e-3,"mA"),(1e-6,"µA"),(1e-9,"nA"),(1e-12,"pA"),(1e-15,"fA")):
        if m>=sc: return sc,un
    return (1e-15,"fA") if m>0 else (1.,"A")
def grid_axes(fig,items,hspace=0.55,wspace=0.42):
    """Axes for the visible panels of a tab: items = [(key, projection or
    None)] in display order -> {key: axes}.  Up to 3 panels share one row,
    4 make a 2 x 2 grid, 5-6 a 2 x 3 grid."""
    n=len(items)
    if n==0: return {}
    nr,nc=(1,n) if n<=3 else ((2,2) if n==4 else (2,3))
    gs=GridSpec(nr,nc,figure=fig,hspace=hspace,wspace=wspace)
    out={}
    for k,(key,proj) in enumerate(items):
        r,c=divmod(k,nc)
        out[key]=fig.add_subplot(gs[r,c],projection=proj) if proj else fig.add_subplot(gs[r,c])
    return out

def row_axes(fig,items,wspace=0.35):
    """One row of panels: items = [(key, width ratio, projection or None)]."""
    if not items: return {}
    gs=GridSpec(1,len(items),figure=fig,wspace=wspace,width_ratios=[w for _,w,_ in items])
    return {key:(fig.add_subplot(gs[k],projection=proj) if proj else fig.add_subplot(gs[k]))
            for k,(key,_,proj) in enumerate(items)}

def sax(ax,xl="",yl="",t="",log=False):
    ax.set_xlabel(xl,fontsize=8); ax.set_ylabel(yl,fontsize=8)
    ax.set_title(t,fontsize=9,fontweight="bold",pad=4); ax.grid(True,lw=0.35,ls="--")
    if log: ax.set_yscale("log")
def sax3(ax,xl="",yl="",zl="",t=""):
    try: ax.locator_params(nbins=4)
    except Exception: pass
    try:                       # a short (thin) box edge gets two ticks, not a pile of them
        from matplotlib.ticker import MaxNLocator
        asp=np.asarray(ax.get_box_aspect(),float); asp=asp/asp.max()
        for k,axis in enumerate((ax.xaxis,ax.yaxis,ax.zaxis)):
            if asp[k]<0.35: axis.set_major_locator(MaxNLocator(2))
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
DEVICE_FILE_VERSION=1
def check_device_json(d):
    """Check a device-file dict before anything of it is applied; raises
    ValueError naming the first bad field (a bad file must leave the editor
    untouched).  Returns a list of warnings (e.g. a newer file version)."""
    def num(v,what,lo=None,hi=None,integer=False,strict_lo=False):
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):
            raise ValueError(f"{what}: expected a number, got {v!r}")
        if integer and int(v)!=v: raise ValueError(f"{what}: expected a whole number, got {v!r}")
        if lo is not None and (v<=lo if strict_lo else v<lo): raise ValueError(f"{what} = {v!r} is too small")
        if hi is not None and v>hi: raise ValueError(f"{what} = {v!r} is too large")
        return v
    if not isinstance(d,dict): raise ValueError("not a SemiSim 3D device file")
    warn=[]
    ver=d.get("version",1)
    if isinstance(ver,(int,float)) and ver>DEVICE_FILE_VERSION:
        warn.append(f"file version {ver} is newer than this program's {DEVICE_FILE_VERSION}: unknown fields are ignored")
    for k in ("Lx","Ly","Lz","Nx","Ny","Nz","regs","cons"):
        if k not in d: raise ValueError(f"missing field '{k}'")
    for k in ("Lx","Ly","Lz"): num(d[k],k,0.,None,strict_lo=True)
    for k in ("Nx","Ny","Nz"): num(d[k],k,2,100000,integer=True)
    for nm in ("regs","cons"):
        if not isinstance(d[nm],list): raise ValueError(f"'{nm}' must be a list")
    for i,r in enumerate(d["regs"]):
        if not isinstance(r,dict): raise ValueError(f"region {i+1}: not an object")
        w=f"region {i+1} ({r.get('label','')})"
        if r.get("mat") not in MATS: raise ValueError(f"{w}: unknown material {r.get('mat')!r}")
        for a in ("x","y","z"):
            v0=num(r.get(a+"0"),f"{w} {a}0"); v1=num(r.get(a+"1"),f"{w} {a}1")
            if not v1>v0: raise ValueError(f"{w}: {a}1 = {v1} must be larger than {a}0 = {v0}")
        if r.get("dtype") not in ("n","p","i",None,""): raise ValueError(f"{w}: doping type must be n, p or i")
        num(r.get("doping",0.),f"{w} doping",0.,1e22)
    for i,c in enumerate(d["cons"]):
        if not isinstance(c,dict): raise ValueError(f"contact {i+1}: not an object")
        w=f"contact {i+1} ({c.get('label','')})"
        num(c.get("face"),f"{w} face",0,len(Con.FNAME)-1,integer=True)
        num(c.get("bc"),f"{w} type (bc)",0,len(Con.BCNAME)-1,integer=True)
        num(c.get("V",0.),f"{w} voltage",-1e5,1e5)
    if len(d["cons"])>32: raise ValueError(f"{len(d['cons'])} contacts: the solver holds at most 32")
    so=d.get("solver") or {}
    if "T" in so: num(so["T"],"temperature (solver T)",20.,2000.)
    ms=(d.get("mesh") or {}).get("type","rect")
    if ms not in ("rect","tri"): raise ValueError(f"unknown mesh type {ms!r} (rect or tri)")
    return warn

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
                 tox=10.,nox=0.,k0p=0.,k1p=1.,net="",stack=None,shape=None,metal=""):
        self.face=face; self.i0p=i0p;self.i1p=i1p;self.j0p=j0p;self.j1p=j1p
        # contact metal (METALS key; "" = ideal ohmic / Schottky with phi_m = pm)
        self.metal=metal if metal in METALS else ""
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

def con_full(face,bc,V,lbl="",pm=4.05,tox=10.,nox=0.,net="",stack=None,metal=""):
    return Con(face,0.,1.,0.,1.,bc,V,pm,lbl,tox,nox,net=net,stack=stack,metal=metal)
def con_part(face,i0p,i1p,j0p,j1p,bc,V,lbl="",pm=4.05,tox=10.,nox=0.,net="",stack=None,metal=""):
    return Con(face,i0p,i1p,j0p,j1p,bc,V,pm,lbl,tox,nox,net=net,stack=stack,metal=metal)
def con_box(x0,y0,z0,x1,y1,z1,L,bc,V,lbl="",pm=4.05,net="",shape=None,metal=""):
    """Buried electrode occupying the box (nm) in a domain of size L=(Lx,Ly,Lz).
    A Gate box is a metal/poly electrode (e.g. inside an oxide-lined trench);
    an Ohmic box on semiconductor is a buried contact.  shape: optional
    cross-section outline cutting the box (see Reg)."""
    Lx,Ly,Lz=L
    return Con(6,x0/Lx,x1/Lx,y0/Ly,y1/Ly,bc,V,pm,lbl,10.,0.,z0/Lz,z1/Lz,net=net,shape=shape,metal=metal)

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
    cons=[con_full(0,1,0.,"Schottky Metal (W)",4.72,metal="W"),con_full(1,0,0.,"Ohmic Back")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-1.5,ivE=0.7,ivN=45,
                dense_x=[0.,450.],dense_y=[],dense_z=[],
                desc="W/n-Si Schottky diode (phi_B 0.67 eV) — dense near metal interface")

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
                ivc=1,ivS=0.,ivE=3.,ivN=31,wl=200./320.,models=255,
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
                ivc=1,ivS=-3.,ivE=0.,ivN=31,wl=200./320.,models=255,
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
    Titanium Schottky anode on the n- drift (Ti/4H-SiC: phi_B 1.10 eV,
    measured) between p+ stripes that pinch off the Schottky regions under
    reverse bias (lower leakage than a plain Schottky diode).  The p+ stripes
    have their own ohmic contacts (in a real part a separate p-ohmic metal
    under the anode); all anode contacts are one net, 'Anode'.  A Schottky
    metal directly on p+ SiC would be a 2.1 eV hole barrier and leave the
    grid floating.

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

    # anode: Ti Schottky segments on the n- surface, ohmic segments on the p+
    edges=[0.]+[e for cx in pp_x_centers for e in (cx-pp_w/2,cx+pp_w/2)]+[Lx]
    cons=[]
    for k in range(0,len(edges),2):
        cons.append(con_part(2,edges[k]/Lx,edges[k+1]/Lx,0.,1.,1,0.,f"Schottky (Ti) {k//2+1}",
                             4.5,net="Anode",metal="Ti"))
    for k,cx in enumerate(pp_x_centers):
        cons.append(con_part(2,(cx-pp_w/2)/Lx,(cx+pp_w/2)/Lx,0.,1.,0,0.,f"p+ ohmic {k+1}",net="Anode"))
    cons.append(con_full(3, 0, 0., "Cathode"))
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-1200.,ivE=2.,ivN=31,bv_thr=0.05,
                dense_x=[c for c in pp_x_centers] + [c-pp_w/2 for c in pp_x_centers]
                       + [c+pp_w/2 for c in pp_x_centers],
                dense_y=[0.,pp_d,sub_top],
                dense_z=[],
                desc="4H-SiC JBS diode — 1.2kV, Ti Schottky (1.1 eV) + ohmic p+ JBS stripes")

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
#   * Impact ionisation (avalanche) is OFF in these templates: switch it on
#     in the Solver page to get a breakdown voltage (a local-field model;
#     with it off the GUI's breakdown marker is a heuristic knee detector).
#     Blocking-state fields and depletion widths are computed either way and
#     can be compared with Ec.
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
    Molybdenum barrier: phi_B 0.68 eV on n-Si (measured), the range of
    40 V rectifier barriers.  Epi: Nd=1.2e16 cm^-3, 3 um -> BV_ideal ~ 46 V.
    The cell is 1 um x 1 um; a 1N5819 die is ~1.5 mm x 1.5 mm (x 2.25e6)."""
    Lx,Ly,Lz = 1000., 4500., 1000.; Nx,Ny,Nz = 8, 120, 6
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,        "n",1.2e16,"n epi (40V)"),
        Reg("Si",0,3000,0,Lx,Ly,Lz,     "n",1e19,  "n+ substrate"),
    ]
    cons=[con_full(2,1,0.,"Schottky Anode (Mo)",4.73,metal="Mo"), con_full(3,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=-50.,ivE=0.6,ivN=41,bv_thr=0.05,
                dense_x=[],dense_y=[],dense_z=[],
                desc="Si Schottky rectifier 40V/1A (1N5819-like) - Mo barrier, phi_B 0.68 eV")

def T_rw_zener_5v1():
    """Si Zener voltage reference, 5.1 V class (1N4733A-like).
    Both sides degenerate (p+ 5e19 / n+ 2e19) so the junction is only ~10 nm
    wide and the field reaches >5 MV/cm at 5 V - the tunnelling regime.
    NOTE: no band-to-band tunnelling model, so the 5.1 V knee is NOT
    reproduced. What IS correct and instructive is the ultra-thin, ultra-high
    field depletion layer that makes Zener action possible.  Leave impact
    ionisation off here: the local avalanche model ignores the dead space of
    a 10 nm junction (carriers cannot gain the ionisation energy in so short
    a distance) and would show an avalanche that a real junction this thin
    does not have."""
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

def T_metal_contacts():
    """Metal contacts: ohmic or rectifying - aluminium on silicon.
    Two aluminium contacts on an n-type layer (2e16 cm^-3, 1.5 um): the left
    one on a shallow n+ implant (1e20, 120 nm), the right one directly on the
    n- layer.  Both are plain 'Ohmic' contacts with Metal = Al: the metal puts
    the measured Al/n-Si Schottky barrier (0.72 eV) at BOTH interfaces, and
    the physics of the interface decides what each contact does (Solver page:
    metal contacts = thermionic emission over the image-force-lowered barrier
    + tunnelling through it):

    - on 1e20 the depletion barrier is only ~3 nm thick: electrons tunnel
      straight through it (field emission) - specific contact resistivity
      ~1e-6 Ohm cm2, an ohmic contact;
    - on 2e16 it is ~0.2 um thick, tunnelling is negligible and only
      thermionic emission over the barrier remains - a Schottky diode.

    Sweep: the right contact from -2 V to +0.5 V, the left one grounded.
    Simulated: forward 4 nA at 0.3 V with an ideality of 1.04 (the
    image-force lowering shrinks with forward bias) until the n- path's
    resistance bends the curve above 0.4 V; reverse -0.06 pA at -0.3 V
    rising to -0.13 pA at -2 V - the 'soft' reverse characteristic of real
    Schottky diodes (image-force lowering and thermionic-field emission grow
    with the field; the ideal-diode model saturates flat).  Rectification
    6e4 at 0.3 V.

    Try:
    - lower the n+ doping to 1e18 (Regions page): the left contact turns
      into a second barrier (back-to-back diodes) and I(+0.3 V) falls ~2000x;
    - Metal = ErSi2 (0.28 eV) for both: a leaky, almost ohmic contact even on
      the n- layer (uA in reverse instead of pA);
    - Metal = Pt (0.90 eV) on the n- contact: 1000x less reverse current;
    - untick the metal-contact physics (Solver page): both contacts become
      ideal ohmic and the device a resistor.
    Currents are for a 200 nm deep slice (Lz)."""
    Lx,Ly,Lz=3000.,1500.,200.; Nx,Ny,Nz=60,36,3
    regs=[Reg("Si",0,0,0,Lx,Ly,Lz,      "n",2e16,"n- Si layer"),
          Reg("Si",0,0,0,900,120,Lz,    "n",1e20,"n+ contact implant")]
    cons=[con_part(2,100/Lx,800/Lx,0.,1.,0,0.,"Al on n+",metal="Al"),
          con_part(2,2000/Lx,2900/Lx,0.,1.,0,0.,"Al on n-",metal="Al")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=-2.,ivE=0.5,ivN=26,
                dense_x=[100.,800.,900.,2000.,2900.],dense_y=[0.,120.],dense_z=[],
                desc="Al contacts on n+ (ohmic, tunnelling) and n- Si (Schottky) - same metal")

def T_defect_diode():
    """Defects: leaky p-n junction - a silicide particle, a void and a dislocation.
    A small 3-D silicon p+/n diode (2 x 2 x 1 um: p+ 1e19 on top, n- 1e16,
    n+ substrate; anode on top, cathode underneath) made from float-zone
    quality silicon (Solver page: background SRH lifetime 1 ms) with three
    common process defects in the plane z = 500 nm (Defects page):

    - a NiSi precipitate (80 nm) in the depletion region - nickel
      contamination; its surface states pin the Fermi level at the NiSi/Si
      barrier (0.65 eV) and recombine like a floating Schottky contact;
    - a void (150 nm) with an unpassivated inner surface (D_it 1e12);
    - a decorated dislocation running through the whole device along y
      (1e6 deep states per cm, Ec - 0.45 eV).

    Sweep: anode -2 V to +0.6 V.  The precipitate dominates everything:
    reverse leakage about -1.4e-14 A, nearly flat with bias (its surface
    generates carriers wherever the depletion region surrounds it), and the
    recombination current at low forward bias: I(0.3 V) 1.4e-11 A against
    8e-15 A for the clean diode, an ideality of ~1.7 instead of 1.0.
    Separately (delete the other rows on the Defects page): at +0.3 V the
    dislocation raises the current 7x, the void 12x, the particle 1900x; at
    +0.5 V 2x, 5x and 94x.  Without the particle the reverse current of this
    small diode is below the solver's resolution (~3e-16 A here; the Results
    table flags it).

    Try:
    - Slicer: quantity 'trap density' (where the states sit) or 'log10|R|'
      (recombination at +0.4 V, generation at -2 V);
    - move the particle into the neutral n-region (y = 1500): its reverse
      leakage disappears (no depletion region around it), while it still
      recombines injected holes in forward bias (3.8e-13 A at 0.3 V);
    - Defects page, Random scatter: copper precipitates (Cu3Si) with a density
      of 1e12 cm^-3 in the n- layer - a contaminated wafer."""
    Lx,Ly,Lz=2000.,2000.,1000.; Nx,Ny,Nz=34,46,16
    regs=[Reg("Si",0,0,0,Lx,Ly,Lz,"n",1e16,"n- Si (1e16)"),
          Reg("Si",0,0,0,Lx,300,Lz,"p",1e19,"p+ anode (1e19)"),
          Reg("Si",0,1800,0,Lx,Ly,Lz,"n",1e19,"n+ substrate")]
    cons=[con_full(2,0,0.,"Anode"),con_full(3,0,0.,"Cathode")]
    defects=[SD.Defect("Metal particle (silicide)",x=600,y=550,z=500,a=80.,metal="NiSi",label="NiSi precipitate"),
             SD.Defect("Void",x=1150,y=900,z=500,a=150.,label="void"),
             SD.Defect("Dislocation (Si)",x=1600,z=500,axis="y",density=1e6,label="dislocation")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,defects=defects,tau0=1e-3,
                ivc=0,ivS=-2.,ivE=0.6,ivN=27,view=("XY",500.),cut="Y",
                dense_x=[],dense_y=[300.,1800.],dense_z=[],
                desc="Si p+/n diode (FZ quality) with a NiSi precipitate, a void and a dislocation - leakage and recombination")

def T_defect_tid():
    """Defects: NMOS after irradiation - total ionizing dose (oxide charge + D_it).
    The NMOS template after a dose of ionizing radiation (Defects page,
    'Ionizing dose'): holes trapped in the gate oxide near the interface
    (+5e11 q/cm2) and radiation-generated interface traps (D_it 5e11
    cm^-2 eV^-1, acceptor-like above midgap, donor-like below).

    Sweep: gate 0 to 1.5 V at V_D = 50 mV (transfer curve; Drain current).
    Pristine (delete the defect): V_th 0.574 V at constant current
    (1e-7 A x W/L; the I-V tab also gives the max-g_m extrapolation, 0.618 V),
    swing 83.9 mV/dec, I_off 1.9e-14 A.  Irradiated: V_th 0.501 V (0.547 V),
    swing 92.9 mV/dec, I_off 3.4e-13 A (18x).  The swing grows by
    kT/q ln10 C_it/C_ox = 9.7 mV/dec, and the threshold shift is a
    difference: the oxide charge alone moves it by -qQ/C_ox = -162 mV, but at
    threshold the acceptor-like traps between midgap and the Fermi level are
    filled (negative) and give back about +90 mV - the classic 'rebound' of
    irradiated NMOS (in a PMOS both act in the same direction).

    Try: twice the dose (D_it 1e12, charge 1e12): V_th 0.430 V, 99.5 mV/dec,
    I_off 220x; only the interface traps (charge 0): V_th 0.664 V (+90 mV),
    the same 92.9 mV/dec; only the oxide charge (D_it 0): -163 mV."""
    t=T_nmos()
    t["cons"][1].V=0.05                  # Drain: linear region
    t["defects"]=[SD.Defect("Ionizing dose (MOS: oxide charge + D_it)",label="ionizing dose")]
    t.update(ivc=2,ivS=0.,ivE=1.5,ivN=31,ivnet="Gate",ivout="Drain",
             desc="NMOS after total ionizing dose: +5e11 q/cm2 oxide charge and D_it 5e11 - V_th shift and swing")
    return t

def T_lateral_pnp():
    """Lateral PNP of a junction-isolated bipolar IC (the classic
    uA741 / LM-series level-shifter device), 2-D cross-section.

    The base is the n-epitaxial layer (5e15, 8 um) of the NPN process; the
    emitter and a collector RING are two p+ base-diffusions (1e19, 2 um deep)
    side by side at the surface, 2 um apart - that gap is the LATERAL base
    width.  The section cuts the ring on both sides of the emitter (contacts
    'Collector L/R' on the net 'Collector').  An n+ buried layer (NBL, 1e19)
    under the device reflects holes injected downwards from the emitter
    bottom (retarding field) back towards the collector and suppresses the
    parasitic VERTICAL PNP emitter -> epi -> p-substrate; p+ isolation walls
    tie to the substrate, which sits at the lowest potential (-5 V).
    n+ base contact (emitter diffusion) outside the ring.

    Bias: base 0 V, collector -2 V, substrate -5 V; the sweep is the emitter
    (V_EB 0.3 -> 0.8 V): a Gummel plot of all terminals.  Simulated, as for
    the textbook device: beta 40-50 at low current, falling to ~30 by
    40 uA (the lightly doped base goes into high injection early), and a
    substrate current of ~2 % of I_C.  Delete the NBL region to watch the
    vertical PNP: the substrate current rises to ~14 % of I_C (beta goes
    UP to 60-75, because holes that recombined in the NBL now escape to the
    substrate instead).  Currents are per 1 um of ring perimeter (Lz)."""
    Lx,Ly,Lz = 24000., 12000., 1000.; Nx,Ny,Nz = 64, 52, 3
    epi, nbl0, nbl1 = 8000., 6500., 9000.
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,            "p",1e15,"p substrate"),
        Reg("Si",0,0,0,Lx,epi,Lz,           "n",5e15,"n-epi (base)"),
        Reg("Si",2000,nbl0,0,22000,nbl1,Lz, "n",1e19,"n+ buried layer"),
        Reg("Si",0,0,0,2000,epi,Lz,         "p",1e18,"p+ isolation"),
        Reg("Si",22000,0,0,Lx,epi,Lz,       "p",1e18,"p+ isolation"),
        Reg("Si",2500,0,0,4000,1000,Lz,     "n",1e20,"n+ base contact"),
        Reg("Si",5000,0,0,8500,2000,Lz,     "p",1e19,"p+ collector (ring)"),
        Reg("Si",10500,0,0,13500,2000,Lz,   "p",1e19,"p+ emitter"),
        Reg("Si",15500,0,0,19000,2000,Lz,   "p",1e19,"p+ collector (ring)"),
    ]
    fx=lambda a,b:(a/Lx,b/Lx)
    cons=[
        con_part(2,*fx(10500,13500),0.,1.,0,0.7,"Emitter"),
        con_part(2,*fx(2500,4000),0.,1.,0,0.,"Base"),
        con_part(2,*fx(5000,8500),0.,1.,0,-2.,"Collector L",net="Collector"),
        con_part(2,*fx(15500,19000),0.,1.,0,-2.,"Collector R",net="Collector"),
        con_full(3,0,-5.,"Substrate"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=0,ivS=0.3,ivE=0.8,ivN=26,ivnet="Emitter",ivout="Collector",cut="X",
                dense_x=[5000.,8500.,10500.,13500.,15500.,19000.],dense_y=[2000.,epi],dense_z=[],
                desc="Lateral PNP (junction-isolated bipolar IC) - p+ emitter + collector ring in n-epi, NBL; Gummel plot")

def T_soi_lateral_npn():
    """SOI lateral NPN - the symmetric lateral bipolar transistor on thin SOI
    (T. H. Ning et al., IBM: CMOS-compatible bipolar in a 60-100 nm film).

    Emitter, base and collector lie side by side in a 60 nm Si film on
    150 nm buried oxide: n+ emitter (1e20) | p base (2e18, W_B = 100 nm) |
    n+ collector (1e20) - symmetric, no buried layer, no vertical parasitic.
    The base is contacted from the top through a p+ polysilicon extrinsic
    base (30 nm thick, 40 nm wide, 1e20, oxide spacers on both sides): its
    p+/p high-low junction keeps the electrons crossing the base away from
    the base metal - an ohmic contact straight on a 100 nm base would
    swallow them (beta < 1).  Emitter and collector contacts on the film
    ends; the p handle wafer under the BOX is grounded.  V_CE = 1 V; the
    sweep is V_BE 0.5 -> 1.0 V (Gummel plot).

    Simulated: an ideal collector current (60 mV/dec) over six decades and
    beta = 16-17 up to ~1 uA, falling at high injection.  About 60 % of
    the base current is electrons that still reach the extrinsic base, the
    rest is hole back-injection into the emitter (band-gap narrowing on):
    the width of the p+ link trades base resistance against gain
    (60 nm -> beta 11, 20 nm -> beta 23 - edit the region to try).
    Currents are per 100 nm of device width (Lz)."""
    Lx,Ly,Lz = 600., 340., 100.; Nx,Ny,Nz = 64, 44, 3
    tpoly, tsi, tbox = 30., 60., 150.
    ys, yb = tpoly, tpoly+tsi                 # film top / bottom
    regs=[
        Reg("SiO2",0,0,0,Lx,ys,Lz,         "i",0.,  "oxide spacers / passivation"),
        Reg("Si",280,0,0,320,ys,Lz,        "p",1e20,"p+ poly extrinsic base"),
        Reg("Si",0,yb+tbox,0,Lx,Ly,Lz,     "p",1e15,"p handle wafer"),
        Reg("SiO2",0,yb,0,Lx,yb+tbox,Lz,   "i",0.,  "buried oxide"),
        Reg("Si",0,ys,0,250,yb,Lz,         "n",1e20,"n+ emitter"),
        Reg("Si",250,ys,0,350,yb,Lz,       "p",2e18,"p base (100 nm)"),
        Reg("Si",350,ys,0,Lx,yb,Lz,        "n",1e20,"n+ collector"),
    ]
    cons=[
        Con(0,ys/Ly,yb/Ly,0.,1.,0,0.,  4.05,"Emitter"),
        Con(1,ys/Ly,yb/Ly,0.,1.,0,1.0, 4.05,"Collector"),
        con_part(2,282./Lx,318./Lx,0.,1.,0,0.8,"Base"),
        con_full(3,0,0.,"Handle"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=2,ivS=0.5,ivE=1.0,ivN=26,ivnet="Base",ivout="Collector",cut="X",
                dense_x=[250.,280.,320.,350.],dense_y=[ys,yb],dense_z=[],
                desc="SOI lateral NPN (IBM-type symmetric lateral bipolar) - 60nm film, 100nm base, p+ poly extrinsic base; Gummel plot")

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

def _comb_bjt(pol):
    """Vertical power BJT section with an interdigitated ('comb') emitter, full
    3-D: x across the fingers, y = depth (top surface y = 0), z along them.
    In x the section runs from the middle of one emitter finger to the
    middle of the next-but-one (half fingers on the mirror planes x = 0 and
    x = Lx, a whole one in the centre, two base-contact fingers between);
    in z from the base bus (z = 0) to the emitter bus (z = Lz), so the
    emitter and the base contact are two combs whose fingers interleave."""
    Lx,Ly,Lz = 280000., 50000., 200000.; Nx,Ny,Nz = 44, 46, 22
    e,b = ("p","n") if pol=="pnp" else ("n","p")    # emitter & collector type, base type
    s = -1. if pol=="pnp" else 1.
    wE,wB,g,bus = 80000., 30000., 15000., 20000.    # emitter / base-contact finger width, gap, bus width
    yBC,yE,yB1,yB2,ySub = 2500., 5000., 8000., 12000., 42000.
    NB1 = 2e16 if pol=="pnp" else 6e16              # active base (NPN: higher - electrons are faster)
    efin=[(0.,wE/2),(Lx/2-wE/2,Lx/2+wE/2),(Lx-wE/2,Lx)]
    bfin=[(wE/2+g,wE/2+g+wB),(Lx/2+wE/2+g,Lx/2+wE/2+g+wB)]
    ze0,zb1 = bus+g, Lz-bus-g                       # emitter fingers start / base fingers end
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,        e,5e14,f"{e}- collector drift (30um)"),
        Reg("Si",0,ySub,0,Lx,Ly,Lz,     e,1e19,f"{e}+ substrate (collector)"),
        Reg("Si",0,0,0,Lx,yB2,Lz,       b,1e16,f"{b} base"),
        Reg("Si",0,0,0,Lx,yB1,Lz,       b,NB1, f"{b} base (active)"),
        Reg("Si",0,0,Lz-bus,Lx,yE,Lz,   e,1e19,f"{e}+ emitter bus"),
    ]
    for k,(x0,x1) in enumerate(efin):
        regs.append(Reg("Si",x0,0,ze0,x1,yE,Lz-bus,e,1e19,f"{e}+ emitter finger {k+1}"))
    regs.append(Reg("Si",0,0,0,Lx,yBC,bus,b,1e19,f"{b}+ base bus"))
    for k,(x0,x1) in enumerate(bfin):
        regs.append(Reg("Si",x0,0,bus,x1,yBC,zb1,b,1e19,f"{b}+ base finger {k+1}"))
    ins=5000.                                       # contact windows 5 um inside the diffusions
    fx=lambda a: a/Lx
    fz=lambda a: a/Lz
    cons=[con_part(2,0.,1.,fz(Lz-bus+ins),1.,0,0.,"Emitter bus",net="Emitter")]
    for k,(x0,x1) in enumerate(efin):
        cons.append(con_part(2,fx(x0+ins if x0>0 else 0.),fx(x1-ins if x1<Lx else Lx),fz(ze0+ins),fz(Lz-bus),
                             0,0.,f"Emitter finger {k+1}",net="Emitter"))
    cons.append(con_part(2,0.,1.,0.,fz(bus-ins),0,0.75*s,"Base bus",net="Base"))
    for k,(x0,x1) in enumerate(bfin):
        cons.append(con_part(2,fx(x0+ins),fx(x1-ins),fz(bus),fz(zb1-ins),0,0.75*s,f"Base finger {k+1}",net="Base"))
    cons.append(con_full(3,0,5.*s,"Collector"))
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=len(cons)-2,ivS=0.45*s,ivE=(0.85 if pol=="pnp" else 0.80)*s,ivN=9 if pol=="pnp" else 8,
                ivnet="Base",ivout="Collector",
                cut="Y",view=("XZ",1000.),
                dense_x=sorted({v for xx in efin+bfin for v in xx if 0<v<Lx}),
                dense_y=[yBC,yE,yB1,yB2,ySub],dense_z=[bus,ze0,zb1,Lz-bus],
                desc=f"Si {pol.upper()} power BJT 250V/16A (MJ2119{3 if pol=='pnp' else 4}-like) - comb (interdigitated) emitter, 3-D section; Gummel plot (heavy: ~30 s per point)")

def T_rw_bjt_mj21193():
    """Si PNP power BJT, 250 V / 16 A / 250 W (MJ21193-like), comb structure.

    An audio-amplifier output transistor: a vertical PNP whose emitter is a
    comb.  The section (0.28 x 0.20 mm, full 3-D) holds one whole p+ emitter
    finger (80 um wide) and two halves, joined by an emitter bus at the far
    end, interleaved with two n+ base-contact fingers (30 um, 15 um gaps)
    that hang from a base bus at the near end.  Every part of the emitter is
    within 55 um of a base contact: the base current flows sideways under the
    emitter, and its voltage drop would otherwise turn off the middle of a
    wide emitter (current crowding at the edges).  onsemi builds the real
    part with a 'Perforated Emitter' instead - one large emitter with an array
    of holes for the base contact - which serves the same purpose.

    Vertical structure: p+ emitter (1e19, 5 um), n base (2e16 to 8 um, 1e16 to
    12 um: ~7 um of base under the emitter, for an f_T of a few MHz), p-
    collector drift 5e14 x 30 um (a punch-through design: the depletion
    reaches the substrate near 370 V, and at V_CB = 400 V the peak field is
    2.4e5 V/cm, just under the avalanche field of such silicon - an ideal
    BV_CBO a little above the 400 V rating), p+ substrate = collector
    contact at the bottom.  Open the Structure tab: the section
    plane starts as a top view (XZ, 1 um deep) that shows the two combs.

    Bias: V_CE = -5 V, V_BE = -0.75 V; the sweep is the base, -0.45 to
    -0.85 V (Gummel plot).  Simulated: hFE 44 and flat from nanoamps to
    about a milliamp per section, 42 at 3.8 mA, 35 at 17 mA and 17 at 39 mA,
    where the lightly doped collector goes into quasi-saturation (its ohmic
    drop forward-biases the collector junction: see the hole density in the
    drift region in the 3D Slicer).  Currents are per section: a 5 x 5 mm
    die - the area the datasheet's C_ob <= 500 pF at V_CB = 10 V allows - is
    ~450 sections, so 17 mA here is ~8 A and 39 mA ~17 A; the datasheet asks
    hFE 25-75 at 8 A (V_CE = 5 V) and >= 8 at 16 A.  Not modelled:
    second breakdown and self-heating (the SOA), AC; impact ionisation
    (V_CEO(sus)) is off - the Solver page switches it on - and f_T is the
    I-V tab's quasi-static estimate.  Heavy: ~20-30 s per bias point on two
    cores (a mesh 3x finer moves the currents by ~6 %, hFE by ~2 %)."""
    return _comb_bjt("pnp")

def T_rw_bjt_mj21194():
    """Si NPN power BJT, 250 V / 16 A / 250 W (MJ21194-like), comb structure.

    The complement of the MJ21193 template - the same comb section, layers
    and contacts with every doping type swapped: n+ emitter comb, p base,
    n- collector drift (5e14 x 30 um), n+ substrate.  The active base is
    doped higher (6e16 instead of 2e16): electrons cross a p base about
    2.5x faster than holes cross an n base, so the same base would give the
    NPN twice the gain - complementary pairs are matched by design like this.

    Bias: V_CE = 5 V, V_BE = 0.75 V; the sweep is the base, 0.45 to 0.80 V
    (Gummel plot).  Simulated: hFE 56 and flat up to ~20 mA per section
    (~10 A for the 5 x 5 mm die, ~450 sections), still 46 at 74 mA (~33 A):
    the n- drift conducts electrons ~3x better than the PNP's p- drift
    conducts holes, so quasi-saturation sets in at higher currents than in
    the PNP.
    Datasheet: hFE 25-75 at 8 A, V_CE = 5 V.  Not modelled: second
    breakdown, self-heating, AC (impact ionisation is off: Solver page).
    Heavy: ~30 s per bias point."""
    return _comb_bjt("npn")

def T_rw_ldo_lateral_pnp():
    """Lateral PNP pass transistor of a 1 A bipolar LDO regulator (LM2940-class), one emitter cell.

    Bipolar low-dropout regulators of the LM2940 class pass the load current
    through a PNP transistor, and in a junction-isolated bipolar process
    that is a large LATERAL PNP: many small round p emitters (0.5 mil =
    12.7 um across), each ringed by a p collector a few um away, in the
    n-epi, with an n+ buried layer under the array.  National
    Semiconductor's patents for such regulators describe 350 emitters, ~3 mA
    each, for 1 A with 0.5 V dropout - the LM2940's ratings.  The template
    is one emitter cell, simulated as its quarter (mirror planes x = 0 and
    z = 0 through the emitter centre), on the triangular-prism mesh - the
    outlines are true circles.  The section plane starts as a top view (XZ,
    1 um deep).

    Layers: p substrate (1e15); n-epi 3e15 x 10 um = the base; n+ buried layer
    (1e19, 8-11 um) - it lowers the base resistance and turns back the holes
    the emitter injects downwards, so few reach the substrate; p emitter
    (radius 6.35 um) and p collector ring (radius 11.35-21 um), 6e18 x 3 um
    (the NPN base diffusion), i.e. a 5 um lateral base; the base is
    contacted at the cell corner through an n+ sinker that reaches the
    buried layer (the NPN's collector reach-through) - with a shallow n+
    contact instead, the lightly doped epi between contact and buried layer
    debiases the emitter by tens of mV.

    Bias - an LDO at full load: base 0 V, emitter +0.85 V (V_in), collector
    2 V below the emitter (V_out), substrate (ground) at -5 V.  Run gives
    0.69 mA per quarter cell = 2.8 mA per emitter at beta 25: 1 A from 350
    emitters takes ~40 mA of base drive, the range of the LM2940's ground
    current at 1 A (30 mA typical, 45 mA maximum).  Beta falls steeply with
    current - 200 at 0.15 mA, 56 at 1.25 mA, 17 at 4.1 mA per emitter: high
    injection in the 3e15 base - which is why these regulators draw tens of
    mA at full load.  The sweep is the dropout: the collector from 2 V to
    0.05 V below the emitter at fixed base drive.  The current holds to
    within 6 % down to V_EC = 0.2 V, then the transistor saturates: the
    collector junction forward-biases, I_C falls by 23 % at 0.1 V and 60 %
    at 0.05 V, while the base current rises x4 and the substrate current
    x10 - the extra ground current every PNP LDO draws in dropout.  The
    datasheet's 0.5 V dropout at 1 A adds the metal, bond-wire and pin
    resistances to the ~0.15 V of the transistor itself.  Currents are per
    quarter cell (x4 per emitter, x1400 for 350 emitters).  The in-plane
    mesh is what matters here (currents within ~1 % and beta within ~4 % of
    a mesh with 2.6x the nodes); Run takes ~40 s on two cores, the sweep
    ~2 min."""
    Lx,Ly,Lz = 30000., 20000., 30000.; Nx,Ny,Nz = 80, 24, 80   # in-plane resolution matters, depth hardly
    Re,WB,Rc2,Rbc = 6350., 5000., 21000., 5000.
    yP,yEpi,yB0,yB1 = 3000., 10000., 8000., 11000.
    VE,VC = 0.85, 0.85-2.0
    emi=SM.Shape.circle("y",0.,0.,Re)
    regs=[
        Reg("Si",0,0,0,Lx,Ly,Lz,      "p",1e15,"p substrate"),
        Reg("Si",0,0,0,Lx,yEpi,Lz,    "n",3e15,"n-epi (base)"),
        Reg("Si",0,yB0,0,Lx,yB1,Lz,   "n",1e19,"n+ buried layer"),
        Reg("Si",0,0,0,Lx,yP,Lz,      "p",6e18,"p collector ring",shape=SM.Shape.circle("y",0.,0.,Rc2)),
        Reg("Si",0,0,0,Lx,yP,Lz,      "n",3e15,"n-epi (lateral base)",shape=SM.Shape.circle("y",0.,0.,Re+WB)),
        Reg("Si",0,0,0,Lx,yP,Lz,      "p",6e18,"p emitter",shape=emi),
        Reg("Si",0,0,0,Lx,yB0,Lz,     "n",1e19,"n+ sinker (base contact)",shape=SM.Shape.circle("y",Lx,Lz,Rbc)),
    ]
    w=0.6*Re                                      # emitter contact window (quarter)
    r0,r1=Re+WB+1500.,Rc2-1500.                   # collector windows: over the ring
    d0,d1=r0/np.sqrt(2)+200.,r1/np.sqrt(2)-200.
    fx=lambda a: a/Lx
    fz=lambda a: a/Lz
    cons=[
        con_part(2,0.,fx(w),0.,fz(w),0,VE,"Emitter"),
        con_part(2,fx(r0),fx(r1),0.,fz(2500.),0,VC,"Collector x",net="Collector"),
        con_part(2,fx(d0),fx(d1),fz(d0),fz(d1),0,VC,"Collector diagonal",net="Collector"),
        con_part(2,0.,fx(2500.),fz(r0),fz(r1),0,VC,"Collector z",net="Collector"),
        con_part(2,fx(Lx-0.6*Rbc),1.,fz(Lz-0.6*Rbc),1.,0,0.,"Base"),
        con_full(3,0,-5.,"Substrate"),
    ]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=Nx,Ny=Ny,Nz=Nz,regs=regs,cons=cons,
                ivc=1,ivS=VE-2.0,ivE=VE-0.05,ivN=21,ivnet="Collector",ivout="Collector",
                mesh="tri",extrude="y",cut="X",view=("XZ",1000.),
                dense_x=[Re,Re+WB,Rc2],dense_y=[yP,yB0,yB1,yEpi],dense_z=[Re,Re+WB,Rc2],
                desc="Lateral PNP pass transistor of a 1A bipolar LDO (LM2940-class) - round emitter cell, quarter; dropout sweep (heavy: ~40 s)")

# ── BJT effects: small 2-D sections that each show one non-ideal effect, with
#    variants for the A/B compare (I-V tab ▸ A/B compare… ▸ Variants) ──────
def _fx_bjt(NC=1e16,NB=1e17,tC=850.,Vbe=0.6,Vce=2.,emitter=True,base=True,ii=False):
    """The planar NPN of the 'NPN BJT' template (x across, y = depth): n+
    emitter at the top left, p+ extrinsic base at the top right, p base
    (0.25 um below the emitter), n collector tC thick on an n+ buried layer
    contacted at the bottom.  emitter = False leaves the emitter out (no
    region, no contact: the collector-base diode); base = False leaves the
    base contact out (the base region floats)."""
    Lx,Lz=2000.,200.; yB=450.; yBL=yB+tC; Ly=yBL+300.
    Ny=70 if tC<=1000 else 90
    regs=[Reg("Si",0,0,0,Lx,Ly,Lz,      "n",NC,  f"n collector ({NC:.0e})"),
          Reg("Si",0,yBL,0,Lx,Ly,Lz,    "n",1e19,"n+ buried layer"),
          Reg("Si",0,0,0,Lx,yB,Lz,      "p",NB,  f"p base ({NB:.0e})")]
    if emitter: regs.append(Reg("Si",0,0,0,900,200,Lz,"n",2e19,"n+ emitter"))
    regs.append(Reg("Si",1400,0,0,Lx,200,Lz,"p",1e19,"p+ extrinsic base"))
    cons=[]
    if emitter: cons.append(con_part(2,0.,900/Lx,0.,1.,0,0.0,"Emitter"))
    if base: cons.append(con_part(2,1400/Lx,1.,0.,1.,0,Vbe,"Base"))
    cons.append(con_full(3,0,Vce,"Collector"))
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=52,Ny=Ny,Nz=4,regs=regs,cons=cons,
                dense_x=[900.,1400.],dense_y=[200.,yB,yBL],dense_z=[],
                models=M_DEFAULT|(M_II if ii else 0),cut="Y")

def _kirk_bjt(NC):
    t=_fx_bjt(NC=NC,Vbe=0.6,Vce=2.)
    t.update(ivS=0.6,ivE=1.1,ivN=26,ivnet="Base",ivout="Collector",
             desc=f"Kirk effect: NPN Gummel plot, collector {NC:.0e} cm^-3 - f_T and beta vs I_C (I-V tab, BJT panel)")
    return t
def T_fx_kirk():
    """BJT effects: Kirk effect - f_T and gain against collector current.

    The planar NPN of the 'NPN BJT' template with a 0.85 um n collector
    (1e16 cm^-3) on its n+ buried layer; V_CE = 2 V, the sweep is the base
    from 0.6 to 1.1 V (a Gummel plot).  The I-V tab's third panel shows the
    current gain beta = I_C/I_B and the quasi-static transit frequency
    f_T = g_m / (2 pi C_in) against I_C (C_in = dQ/dV_BE: the change of the
    holes stored in the device).  At low current f_T rises with I_C (the
    junction capacitances charge faster); it peaks at 6.3 GHz near I_C =
    1.2e-5 A per cell (~6.7e3 A/cm2 over the emitter contact) and then falls:
    the injected electrons outnumber the collector doping, the field at the
    base-collector junction collapses and the base widens into the collector
    (Kirk effect: base push-out at J_K ~ q v_sat N_C = 1.6e4 A/cm2 for 1e16),
    together with high injection in the base.  Beta falls from 36 to 2 over
    the same range.  See the holes spread into the collector in the 3D
    Slicer (log10(p)) at the last sweep point.

    Variants (A/B compare): the same transistor with a 1e17 collector - its
    f_T peaks at 11.5 GHz near 1.6e-5 A (9e3 A/cm2): a heavier collector
    pushes the Kirk effect to higher current (at the price of a lower
    breakdown voltage and a larger collector capacitance).  Capture the two
    with 'Variants → A / → B', then Run A and B.  ~25-40 s per sweep."""
    return _kirk_bjt(1e16)

def _crowd_bjt(nf):
    """Emitter-crowding section (x across, y depth, 2-D): total emitter
    width 20 um as one stripe (nf = 1) or nf narrow fingers, with p+ base
    contacts beside every stripe; same width, layers and doping."""
    Lx,Ly,Lz=43000.,2000.,1000.
    yE,yB,yX,yC=150.,350.,300.,1500.
    if nf==1:
        bases=[(0.,10500.),(32500.,Lx)]; ems=[(11500.,31500.)]
    else:
        w=20000./nf; bases=[]; ems=[]; x=0.
        bw=(Lx-20000.-2*nf*1000.)/(nf+1)
        for k in range(nf):
            bases.append((x,x+bw)); x+=bw+1000.; ems.append((x,x+w)); x+=w+1000.
        bases.append((x,Lx))
    regs=[Reg("Si",0,0,0,Lx,Ly,Lz,"n",3e16,"n collector (3e16)"),
          Reg("Si",0,yC,0,Lx,Ly,Lz,"n",1e19,"n+ substrate"),
          Reg("Si",0,0,0,Lx,yB,Lz,"p",1e17,"p base (1e17, ~10 kOhm/sq under the emitter)")]
    for k,(a,b) in enumerate(ems): regs.append(Reg("Si",a,0,0,b,yE,Lz,"n",2e19,f"n+ emitter {k+1}"))
    for k,(a,b) in enumerate(bases): regs.append(Reg("Si",a,0,0,b,yX,Lz,"p",1e19,f"p+ base {k+1}"))
    cons=[con_part(2,a/Lx,b/Lx,0.,1.,0,0.0,f"Emitter {k+1}",net="Emitter") for k,(a,b) in enumerate(ems)]
    cons+=[con_part(2,(a+500.)/Lx if a>0 else 0.,(b-500.)/Lx if b<Lx else 1.,0.,1.,0,0.6,f"Base {k+1}",net="Base")
           for k,(a,b) in enumerate(bases)]
    cons.append(con_full(3,0,2.0,"Collector"))
    edges=sorted({v for ab in ems+bases for v in ab if 0<v<Lx})
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=150,Ny=48,Nz=2,regs=regs,cons=cons,
                ivS=0.6,ivE=0.98,ivN=20,ivnet="Base",ivout="Collector",prof="Emitter",
                dense_x=edges,dense_y=[yE,yX,yB,yC],dense_z=[],cut="X",models=M_DEFAULT,
                desc=("Emitter current crowding: one 20 um emitter" if nf==1 else
                      f"Emitter current crowding: {nf} x {20/nf:g} um emitter fingers")+
                     " - Gummel plot; Current tab ▸ Contact = Emitter")
def T_fx_crowding():
    """BJT effects: emitter current crowding - one wide emitter vs narrow fingers.

    A 2-D NPN section 43 um wide with 20 um of emitter, as ONE stripe with a
    p+ base contact on each side.  The base current reaches the middle of
    the emitter sideways through the thin p base under it (0.2 um, 1e17:
    ~10 kOhm/sq, a 'pinched' base); its voltage drop lowers V_BE in the
    middle, so at high current the emitter injects mainly at its edges -
    emitter current crowding.  The sweep is the base (0.6 to 0.98 V,
    V_CE = 2 V).  Open the Current tab with 'Contact' = Emitter (set by the
    template): the 'Current under a contact' panel shows the current
    density across the emitter - flat at low current, rising towards the
    edges once the lateral base drop reaches a few kT/q (J_B ~ 5 A/cm2).
    At the last point (0.98 V) the middle of the stripe carries 1.4 % of the
    current density at its edges (J max / mean 12).

    Variants (A/B compare): the same 20 um of emitter as four 5 um fingers
    with a base contact between every two (the base drop falls with the
    square of the width).  Run both: the Gummel plots part at high current -
    at 0.98 V the fingers deliver 3.5x the collector current of the single
    stripe (1.3x at 0.70 V, the same at 0.6 V), because the middle of the
    wide emitter is turned off; each finger crowds far less (centre / edge
    0.17 at 0.98 V, J max / mean 3.5).  This is why power and RF transistors
    use comb (interdigitated) or perforated emitters (see the MJ21193 /
    MJ21194 templates).  ~1 min per sweep."""
    return _crowd_bjt(1)

def _early_bjt(NB):
    t=_fx_bjt(NC=1e16,NB=NB,Vbe=0.75,Vce=0.)
    t.update(ivS=0.,ivE=5.,ivN=26,ivnet="Collector",ivout="Collector",
             desc=f"Early effect: NPN output curve I_C(V_CE) at V_BE = 0.75 V, base {NB:.0e} cm^-3")
    return t
def T_fx_early():
    """BJT effects: Early effect - output curve and Early voltage.

    The planar NPN with a light p base (5e16 cm^-3, 0.25 um below the
    emitter) at a fixed V_BE = 0.75 V; the sweep is the collector, 0 to 5 V
    (output characteristic).  Past saturation (V_CE < ~0.3 V) the collector
    current keeps rising: the base-collector depletion region widens with
    V_CB and eats into the neutral base, so the electrons cross a thinner
    base (base-width modulation).  The I-V tab fits the flat part of the
    curve and reports the Early voltage V_A = I_C / g_o - V_CE (the dotted
    line meets I = 0 at V = -V_A): about 6 V here, with beta ~ 90 at
    V_CE = 2 V - I_C grows by 85 % from 0.2 to 5 V.

    Variants (A/B compare): base 5e16 vs 2e17 cm^-3 - the heavier base is
    barely depleted, so its V_A is ~42 V (I_C +13 % over the same range), but
    its beta falls to ~17 and its I_C to a sixth: fewer electrons get across
    the heavier base.  Every BJT design trades gain against Early voltage
    this way; the HBT template shows the way out.  ~25 s per sweep."""
    return _early_bjt(5e16)

def _aval_diode(ii):
    Lx,Lz=200.,200.; Wp,Wn,Wnp=1000.,4000.,500.; Ly=Wp+Wn+Wnp
    regs=[Reg("Si",0,0,0,Lx,Ly,Lz,"n",1e16,"n (1e16)"),
          Reg("Si",0,0,0,Lx,Wp,Lz,"p",1e19,"p+ anode"),
          Reg("Si",0,Wp+Wn,0,Lx,Ly,Lz,"n",1e19,"n+ cathode")]
    cons=[con_full(2,0,0.,"Anode"),con_full(3,0,0.,"Cathode")]
    return dict(Lx=Lx,Ly=Ly,Lz=Lz,Nx=3,Ny=260,Nz=3,regs=regs,cons=cons,
                ivc=0,ivS=0.,ivE=-56.,ivN=57,bv_thr=0.,cut="Y",
                dense_x=[],dense_y=[Wp,Wp+Wn],dense_z=[],models=M_DEFAULT|(M_II if ii else 0),
                desc="Avalanche breakdown: Si p+/n (1e16) diode, reverse sweep to -56 V"+
                     ("" if ii else " - impact ionisation OFF"))
def T_fx_avalanche_diode():
    """Avalanche breakdown: Si p+/n junction (BV ~ 52 V).

    An abrupt p+ (1e19) / n (1e16, 4 um) / n+ diode, 1-D, swept in reverse
    to -56 V with impact ionisation on (Solver page: Chynoweth law with the
    van Overstraeten - de Man coefficients).  The leakage current is flat
    (below the current resolution, grey) until the multiplication takes off
    a few volts below breakdown; there the solver stops converging, and the
    I-V tab reports BV where 1/I, extrapolated from the last converged
    points, reaches zero: ~52 V.  Past breakdown the current jumps to a
    branch held back only by the series resistance and the space charge of
    the avalanche carriers.  Theory: the
    ionisation integral of the simulated field reaches 1 at 51.9 V, and
    Baliga's closed form for silicon gives 5.34e13 N^-3/4 = 53.4 V.  The
    3D Slicer quantity 'log10 G avalanche' shows where the pairs are made
    (the high-field part of the depletion region).

    Variants (A/B compare): the same diode with impact ionisation off - its
    current never runs away (no breakdown at any voltage): the avalanche is
    the model, not the field.  See validate.py section 17 for the checks
    behind these numbers."""
    return _aval_diode(True)

def _bv_bjt(mode):
    t=_fx_bjt(NC=1e16,NB=1e17,tC=3000.,Vbe=0.,Vce=0.,ii=True,emitter=(mode!="cbo"))
    t.update(Nx=40,Nz=3,ivS=0.,ivE=60.,ivN=31,ivnet="Collector",ivout="Collector",bv_thr=0.,
             ivfloat="Base" if mode=="ceo" else None,
             desc={"cbo":"BV_CBO: the collector-base junction alone (emitter left out)",
                   "ces":"BV_CES: NPN collector sweep, base shorted to the emitter",
                   "ceo":"BV_CEO: NPN collector sweep, base open (floating, I = 0)"}[mode]+" - impact ionisation on")
    return t
def T_fx_bv_bjt():
    """Avalanche breakdown: NPN BV_CBO, BV_CES and BV_CEO.

    The planar NPN with a thick (3 um, 1e16 cm^-3) collector, impact
    ionisation on, the collector swept from 0 to 60 V in 2 V steps.  Three
    variants (A/B compare: pick two) differ only in what the base and
    emitter do:
    - BV_CES (this one): the base shorted to the emitter, both at 0 V.  The
      holes the avalanche makes flow out to the base contact - sideways
      under the emitter, through the thin base - and the voltage they drop
      there turns the emitter edge on: the transistor amplifies the
      multiplication, and the current runs away at ~34 V.
    - BV_CEO: the base open (a floating terminal, I = 0 - the Sweep page's
      'Floating terminal').  All the holes go into the base, forward-bias
      the emitter and return multiplied by the transistor's gain: the
      current runs away once M = 1 + 1/beta, at ~21 V (between 20 and 22 V,
      the base floating up to 0.9 V) - close to Miller's estimate
      BV_CBO / beta^(1/n) with beta ~ 40, n ~ 4.
    - BV_CBO: the collector-base junction on its own (the emitter left out,
      as an open emitter carries no current): ~56 V - a little above the
      52 V of the p+/n diode template, as the 1e17 base takes a share of
      the voltage.
    BV_CBO > BV_CES > BV_CEO: the order of every bipolar transistor's
    voltage ratings.  Below breakdown the leakage of this small section is
    under the current resolution (grey points); the open base's voltage is
    drawn on the semi-log panel's right axis.  ~0.5-2 min per variant, most
    of it in the failed points at breakdown."""
    return _bv_bjt("ces")

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
    hence essentially zero reverse-recovery charge).  The anode is one net:
    Schottky metal (phi_m 5.0 eV, phi_B 1.2 eV) on the n- surface between the
    grid stripes, ohmic contacts on the p+ stripes."""
    Lx,Ly,Lz = 4000., 14500., 500.; Nx,Ny,Nz = 40, 130, 4
    sub_top, pp_d, pp_w = 12000., 500., 500.
    regs=[Reg("4H-SiC",0,0,0,Lx,Ly,Lz,          "n",1.2e16,"n drift (1200V)"),
          Reg("4H-SiC",0,sub_top,0,Lx,Ly,Lz,    "n",5e18,  "n+ substrate")]
    cxs=(700., 2000., 3300.)
    for cx in cxs:
        regs.append(Reg("4H-SiC",cx-pp_w/2,0,0,cx+pp_w/2,pp_d,Lz,
                        "p",1e18,f"p+ JBS grid @{cx:.0f}nm"))
    # anode net: Schottky metal (phi_m 5.0 eV -> phi_B 1.2 eV) on the n-
    # surface, ohmic contacts on the p+ grid
    edges=[0.]+[e for cx in cxs for e in (cx-pp_w/2,cx+pp_w/2)]+[Lx]
    cons=[con_part(2,edges[k]/Lx,edges[k+1]/Lx,0.,1.,1,0.,f"Schottky Anode {k//2+1}",5.00,net="Anode")
          for k in range(0,len(edges),2)]
    cons+=[con_part(2,(cx-pp_w/2)/Lx,(cx+pp_w/2)/Lx,0.,1.,0,0.,f"p+ ohmic {k+1}",net="Anode")
           for k,cx in enumerate(cxs)]
    cons.append(con_full(3,0,0.,"Cathode"))
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
    spread it with gate/source field plates (impact ionisation is off here;
    the Solver page switches it on)."""
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
                ivc=1,ivS=0.,ivE=VDD,ivN=26,ivnet="Vin",ivfloat="Out",ivout="VDD",models=255,
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
                ivc=2,ivS=0.,ivE=0.8,ivN=17,ivnet="Gate",ivout="Drain",models=255,
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
                ivc=2,ivS=ivS,ivE=ivE,ivN=ivN,ivnet="Gate",ivout="Gate",models=255,
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
                ivc=2,ivS=0.,ivE=0.8,ivN=17,ivnet="Gate",ivout="Drain",models=255,
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
                ivc=2,ivS=0.,ivE=0.7,ivN=15,ivnet="Gate",ivout="Drain",models=255,
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
    "Metal contacts: ohmic vs Schottky (Al)":T_metal_contacts,
    "Defects: leaky p-n junction (particle, void, dislocation)":T_defect_diode,
    "Defects: NMOS after irradiation (oxide charge + D_it)":T_defect_tid,
    "NPN BJT":T_npn,"PNP BJT":T_pnp,"NMOS":T_nmos,"PMOS":T_pmos,
    "Lateral PNP (bipolar IC)":T_lateral_pnp,"SOI lateral NPN":T_soi_lateral_npn,
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
    "RW: Si PNP power BJT 250V (MJ21193)": T_rw_bjt_mj21193,
    "RW: Si NPN power BJT 250V (MJ21194)": T_rw_bjt_mj21194,
    "RW: Si lateral PNP LDO 1A (LM2940)":  T_rw_ldo_lateral_pnp,
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
    # ── BJT effects & breakdown (each with A/B variants) ────────
    "BJT effects: Kirk effect (f_T vs I_C)":       T_fx_kirk,
    "BJT effects: emitter current crowding":       T_fx_crowding,
    "BJT effects: Early effect (output curve)":    T_fx_early,
    "Avalanche: Si p+n diode breakdown":           T_fx_avalanche_diode,
    "Avalanche: NPN BV_CBO, BV_CES, BV_CEO":       T_fx_bv_bjt,
}
# Template variants for the A/B compare (I-V tab): label -> template dict
TEMPLATE_VARIANTS = {
    "BJT effects: Kirk effect (f_T vs I_C)": [
        ("collector 1e16 (early Kirk)",lambda:_kirk_bjt(1e16)),
        ("collector 1e17 (late Kirk)", lambda:_kirk_bjt(1e17))],
    "BJT effects: emitter current crowding": [
        ("one 20 um emitter (crowding)",lambda:_crowd_bjt(1)),
        ("four 5 um fingers",          lambda:_crowd_bjt(4))],
    "BJT effects: Early effect (output curve)": [
        ("base 5e16 (low V_A)", lambda:_early_bjt(5e16)),
        ("base 2e17 (high V_A)",lambda:_early_bjt(2e17))],
    "Avalanche: Si p+n diode breakdown": [
        ("impact ionisation on", lambda:_aval_diode(True)),
        ("impact ionisation off",lambda:_aval_diode(False))],
    "Avalanche: NPN BV_CBO, BV_CES, BV_CEO": [
        ("BV_CES (base shorted to the emitter)",lambda:_bv_bjt("ces")),
        ("BV_CEO (base open)",lambda:_bv_bjt("ceo")),
        ("BV_CBO (emitter left out)",lambda:_bv_bjt("cbo"))],
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
    ("Bipolar & thyristors", ["NPN BJT","PNP BJT","Lateral PNP (bipolar IC)","SOI lateral NPN",
        "GaAs HBT","RW: Si NPN BJT 45V (BC547)",
        "RW: Si PNP power BJT 250V (MJ21193)","RW: Si NPN power BJT 250V (MJ21194)",
        "RW: Si lateral PNP LDO 1A (LM2940)",
        "RW: Si thyristor 800V (BT151)"]),
    ("GaN HEMTs", ["GaN HEMT","RW: GaN RF HEMT 28V (CGH40010)","RW: GaN HEMT 650V (GS66508)"]),
    ("Contacts & defects", ["Metal contacts: ohmic vs Schottky (Al)",
        "Defects: leaky p-n junction (particle, void, dislocation)",
        "Defects: NMOS after irradiation (oxide charge + D_it)"]),
    ("Optoelectronics & solar", ["GaAs LED","Si Solar","RW: Si PIN photodiode (BPW34)",
        "RW: Ge photodiode 1550nm","RW: AlGaAs IR LED 870nm","RW: Si PERC solar cell"]),
    ("BJT effects & breakdown (A/B variants)", ["BJT effects: Kirk effect (f_T vs I_C)",
        "BJT effects: emitter current crowding","BJT effects: Early effect (output curve)",
        "Avalanche: Si p+n diode breakdown","Avalanche: NPN BV_CBO, BV_CES, BV_CEO"]),
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
        fx,fy,fz=SD.refine_planes(t.get("defects") or [],(Lx,Ly,Lz))
        dx=sorted(set(list(t.get("dense_x",[]))+ax+fx))
        dy=sorted(set(list(t.get("dense_y",[]))+ay+fy))
        dz=sorted(set(list(t.get("dense_z",[]))+az+fz))
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
        # metal-semiconductor faces (Schottky or metal contacts): the barrier
        # is only a few nm thick on a heavily doped contact, and tunnelling
        # through it needs the band profile there - 0.25 nm first step
        for c in t.get("cons",[]):
            if c.face<=5 and c.bc!=2 and (c.bc==1 or getattr(c,"metal","")):
                ax=c.face//2; pos=0. if c.face%2==0 else Ls[ax]*1e-9
                g[ax]=insert_layer(g[ax],pos,0.25e-9,1.3,+1 if c.face%2==0 else -1)
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
            arr=np.asarray(g[ax],dtype=float); a_lo,a_hi=arr[0],arr[-1]
            for q in (pos_m-0.5e-9,pos_m):
                # (domain ends taken once: arr is unsorted after the first append,
                # which used to drop the second node of every pin)
                if not (a_lo<q<a_hi): continue
                # a pinned node replaces generated nodes closer than 0.05 nm
                # (a graded node 0.01 nm from a pin only degrades conditioning)
                near=(np.abs(arr-q)<0.05e-9)&(np.abs(arr-q)>1e-13)
                near&=~np.isin(arr,list(pinned[ax]))&(arr>a_lo)&(arr<a_hi)
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
        # box-shaped inclusion defects (voids, particles): pinned faces like
        # insulator regions, so the box is exactly its size (half-open rule)
        for ax,vals in enumerate(SD.pin_planes(t.get("defects") or [],tuple(Ls))):
            for v in vals:
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
    # lines added by adaptive mesh refinement (amr_run): nm per axis
    if t.get("amr"):
        g=[np.asarray(xs,dtype=float),np.asarray(ys,dtype=float),np.asarray(zs,dtype=float)]
        for axn,vals in t["amr"].items():
            ax={"x":0,"y":1,"z":2}[axn]; arr=g[ax]
            add=np.array([q*1e-9 for q in vals if arr[0]<q*1e-9<arr[-1]])
            if not len(add): continue
            d=np.min(np.abs(arr[:,None]-add[None,:]),axis=0)
            out=np.sort(np.concatenate([arr,add[d>0.02e-9]]))
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

MAX_MATERIALS, MAX_CONTACTS = 16, 32        # NM_MAX, NC_MAX of the C core
def _addmat(mi,T):
    """Register a material with the core; a clear error past its limit
    (the core would return -1 and every region of it would stay unpainted)."""
    if not (np.isfinite(T) and T>0): raise ValueError(f"temperature must be above 0 K (got {T})")
    r=api_addmat(mi,T)
    if r<0: raise ValueError(f"too many different materials: the solver core holds at most {MAX_MATERIALS} "
                             "(gate-stack and defect materials count too)")
    return r
def _addcon(*a):
    """Register a contact; a clear error past the core's limit (it would drop
    the contact and shift every later contact index)."""
    r=api_addconx(*a)
    if r<0: raise ValueError(f"too many contacts: the solver core holds at most {MAX_CONTACTS}")
    return r

def _gate_stack(cid, con, mids, T):
    """Layer stack of a Robin gate for tunnelling (ignored unless tunnelling is on)."""
    lay=con.layers(); ids=[]
    for m,_ in lay:
        mi=MATS.get(m,8)
        if mi not in mids: mids[mi]=_addmat(mi,T)
        ids.append(mids[mi])
    if ids and min(ids)>=0:
        api_gstack(cid,len(lay),(I*len(lay))(*ids),(D*len(lay))(*[tt*1e-9 for _,tt in lay]))

_PRISM=None          # the prism mesh of the device in the core (None: tensor mesh)

SEMIS = set(m for m in MATS if m not in INSULATORS)
def _eg_of(name): return SEMI_EG.get(name,1.12)
def _barrier_of(metal,semi):
    if metal in METALS and semi in SEMI_CHI: return schottky_barrier(metal,semi)[0]
    return 0.5*_eg_of(semi)

def _node_names(regs, own):
    """Material name of every node from its owner region index (-1: the first region)."""
    names=np.array([r.mat for r in regs]+[regs[0].mat if regs else "Si"],dtype=object)
    o=np.asarray(own).ravel(); return names[np.where(o<0,len(regs),o)]

def _face_gates_tensor(t, xs, ys, zs):
    """[(nodes, areas m^2)] of the semiconductor face nodes under every face
    (oxide) gate on a tensor mesh - interfaces for interface traps/charge."""
    out=[]; n=(len(xs),len(ys),len(zs)); c=(np.asarray(xs),np.asarray(ys),np.asarray(zs))
    fv=[SD.MeshView._fv(a) for a in c]; L=(t["Lx"],t["Ly"],t["Lz"])
    for con in t["cons"]:
        if con.bc!=2 or con.face>5: continue
        ax=con.face//2; e=con_extent(con,*L); idx=[None,None,None]; idx[ax]=np.array([0 if con.face%2==0 else n[ax]-1])
        for q in range(3):
            if q==ax: continue
            a0,a1=span_idx(c[q],e["xyz"[q]][0]*1e-9,e["xyz"[q]][1]*1e-9); idx[q]=np.arange(a0,a1+1)
        I,J,K=np.meshgrid(idx[0],idx[1],idx[2],indexing="ij")
        o=[q for q in range(3) if q!=ax]; G=(I,J,K)
        area=fv[o[0]][G[o[0]]]*fv[o[1]][G[o[1]]]
        out.append(((I+n[0]*(J+n[1]*K)).ravel(),area.ravel()))
    return out

def defect_hints(t):
    """Bounding boxes (nm) of the defects the prism mesher should refine around."""
    L=(t["Lx"],t["Ly"],t["Lz"]); out=[]
    for d in (t.get("defects") or [])[:SD.REFINE_MAX]:
        if d.preset()["geom"] in ("inclusion","line","plane") and d.shape!="everywhere":
            lo,hi=d.extent(L); out.append((tuple(lo),tuple(hi)))
    return out

def _apply_defects(t, view, mat, Nd, Na, mids, T, gate_faces):
    """Map t['defects'] onto the node arrays (material overrides in place)
    and return the SD.DefectData to load after the contacts."""
    L=(t["Lx"],t["Ly"],t["Lz"])
    dd=SD.apply(t["defects"],view,L,SEMIS,INSULATORS,_eg_of,_barrier_of,gate_faces)
    for k,mn in dd.mat.items():
        mi=MATS[mn]
        if mi not in mids: mids[mi]=_addmat(mi,T)
        mat[k]=mids[mi]; Nd[k]=0.; Na[k]=0.
    return dd

def _load_defects(dd):
    """Deep levels and fixed charge of a mapped defect set into the core."""
    if dd is None: return
    TY={"A":-1,"D":1,"N":0}; RF={"Ec":0,"Ev":1,"Ei":2}
    for key,nodes,dens in dd.traps:
        ty,rf,E,sn,sp=key
        sid=api_trapsp(TY[ty],RF[rf],float(E),float(sn),float(sp))
        if sid is None or sid<0: continue
        api_trapadd(sid,len(nodes),_ip(nodes.astype(np.int32)),_dp(dens))
    qn,qv=dd.charge()
    if len(qn): api_qadd(len(qn),_ip(qn.astype(np.int32)),_dp(qv))
    global _DEFECT_NOTES
    _DEFECT_NOTES=list(dd.notes)
_DEFECT_NOTES=[]

def _apply_tau0(t, mids):
    """Background SRH lifetime t['tau0'] (s) for every semiconductor: None =
    the material defaults, 0 = none (only defects' deep levels recombine)."""
    tau0=t.get("tau0")
    if tau0 is None: return
    names={v:k for k,v in MATS.items()}
    for mid,ci in mids.items():
        if names.get(mid) in SEMIS: api_settau(ci,float(tau0),float(tau0))

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
        if m not in mids: mids[m]=_addmat(m,T)
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
    dd=None
    if any(getattr(r,"shape",None) is not None for r in t["regs"]) or t.get("defects"):
        # shaped regions or defects: node by node (a staircase of the outline)
        own=owner_grid(t["regs"],xs,ys,zs)
        mat,Nd,Na=_node_arrays(t["regs"],own,mids)
        if t.get("defects"):
            view=SD.MeshView(_node_names(t["regs"],own),tensor=(xs,ys,zs))
            dd=_apply_defects(t,view,mat,Nd,Na,mids,T,_face_gates_tensor(t,xs,ys,zs))
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
    P3=None; L3=(t["Lx"],t["Ly"],t["Lz"])
    for con in t["cons"]:
        i0,i1,j0,j1,k0,k1=con.to_idx6(xs,ys,zs)
        V=volt.get(con.label,con.V) if volt else con.V
        robin=con.bc==2 and con.face<=5           # oxide (Robin) gate on a face
        tox=con.eot() if robin else con.tox
        bc_c,pm_c,_=contact_core(con,t["regs"],L3,models)
        cid=_addcon(con.face,i0,i1,j0,j1,k0,k1,bc_c,V,pm_c,tox*1e-9,con.nox*1e4)
        if robin and cid>=0: _gate_stack(cid,con,mids,T)
        if cid>=0 and getattr(con,"shape",None) is not None:
            # shaped box contact: its nodes as a list (staircase of the outline)
            if P3 is None: P3=np.meshgrid(np.asarray(zs)*1e9,np.asarray(ys)*1e9,np.asarray(xs)*1e9,indexing="ij")
            m=SM.box_contact_mask(con,P3[2],P3[1],P3[0],L3).ravel()
            nodes=np.nonzero(m)[0]
            if len(nodes): api_cnodes(cid,len(nodes),_ip(nodes),_dp(np.zeros(len(nodes))))
    _load_defects(dd); _apply_tau0(t,mids)
    return xs,ys,zs

def _setup_tri(t, solver, models, T, volt):
    """setup_device on the triangular-prism mesh."""
    global _PRISM
    _PRISM=None
    api_init()
    mids={}
    for r in t["regs"]:
        m=MATS[r.mat]
        if m not in mids: mids[m]=_addmat(m,T)
    xs,ys,zs=make_grid(t,True)
    tm=dict(t); tm["mesh_hints"]=defect_hints(t)
    pm=SM.PrismMesh(tm,(xs,ys,zs),prism_axis(t,t.get("prism")),INSULATORS,thin=t.get("tri_thin",True))
    api_solver(*solver)
    N,E,rp,cj,h,ar,vol,pos=pm.graph()
    if not api_graph(N,E,_ip(rp),_ip(cj),_dp(h),_dp(ar),_dp(vol),_dp(pos)):
        raise MemoryError(f"Could not set up the prism mesh ({N:,} nodes, {E:,} edges).")
    api_models(models)
    own=pm.node_regions()
    mat,Nd,Na=_node_arrays(t["regs"],own,mids)
    dd=None
    if t.get("defects"):
        gf=[pm.contact_nodes(c) for c in t["cons"] if c.bc==2 and c.face<=5]
        view=SD.MeshView(_node_names(t["regs"],own),graph=(N,rp,cj,ar,vol,pos))
        dd=_apply_defects(t,view,mat,Nd,Na,mids,T,gf)
    api_ndata(_ip(mat),_dp(Nd),_dp(Na),_dp(pm.sheet_charge()))
    ins_id={mids[MATS[m]] for m in INSULATORS if MATS[m] in mids}
    kind=np.where(np.isin(mat,list(ins_id)),1,0).astype(np.int8)
    api_clrcon(); robin=[]; L3=(t["Lx"],t["Ly"],t["Lz"])
    for con in t["cons"]:
        V=volt.get(con.label,con.V) if volt else con.V
        rob=con.bc==2 and con.face<=5
        tox=con.eot() if rob else con.tox
        bc_c,pm_c,_=contact_core(con,t["regs"],L3,models)
        cid=_addcon(con.face,0,0,0,0,0,0,bc_c,V,pm_c,tox*1e-9,con.nox*1e4)
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
    _load_defects(dd); _apply_tau0(t,mids)
    _PRISM=pm
    return xs,ys,zs

# ══════════════════════════════════════════════════════════════
#  PULL SOLUTION
# ══════════════════════════════════════════════════════════════
_PULL_Q=[("phi",0,1.),("n",1,1e-6),("p",2,1e-6),("Ec",3,1.),("Ev",4,1.),
         ("Efn",5,1.),("Efp",6,1.),("Jx",7,1e-4),("Jy",8,1e-4),("Jz",9,1e-4),
         ("mn",10,1e4),("mp",11,1e4),("R",12,1e-6),("kind",13,1.),("Nnet",14,1e-6),
         ("Jg",27,1e-4),("Lqn",28,1.),("Lqp",29,1.),
         ("Qtrap",31,1e-6),("Ntrap",32,1e-6),
         ("Gii",33,1e-6)]  # J in A/cm2, quantum potentials in eV, traps in cm^-3, avalanche generation in cm^-3 s^-1
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
    d["Iraw"]=[api_conIraw(c) for c in range(nc)]     # each contact's own flux (no Kirchhoff substitution)
    d["conv"]=bool(api_conv())                         # convergence of THIS solution
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
# defects: deep-level species, per-node trap densities, fixed charge
api_trapsp  = _opt("api3_trap_species",I,I,I,D,D,D,default=-1)             # type,ref,E(eV),sn,sp(m^2)
api_trapadd = _opt("api3_trap_add",I,I,I,IP,DP,default=0)                  # sid,n,nodes,dens(m^-3)
api_qadd    = _opt("api3_charge_add",I,I,IP,DP,default=0)                  # n,nodes,rho(q/m^3)
api_settau  = _opt("api3_set_mat_tau",I,I,D,D,default=0)                    # mat idx, tau_n0, tau_p0 (s)
api_csums   = _opt("api3_carrier_sums",V,DP)                                # out[2]: sum n dV, sum p dV
QE=1.602176634e-19                                                          # elementary charge (C)

def carrier_charge():
    """Electron and hole charge stored in the semiconductor, (q sum n dV,
    q sum p dV) in C per simulated cell - NaN with a core that predates it."""
    b=(D*2)(float("nan"),float("nan")); api_csums(b)
    return b[0]*QE,b[1]*QE
# model switches (api3_set_models): bit 0..5
M_BGN,M_TAU,M_VSAT,M_QC,M_LB,M_TUN = 1,2,4,8,16,32
M_TE,M_SBT = 64,128        # metal contacts: thermionic emission + image force; barrier tunnelling
M_II = 256                 # impact ionisation (avalanche generation, Chynoweth law)
M_DEFAULT = M_BGN|M_TAU|M_VSAT|M_TE|M_SBT

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

def solve_floating(fidx,V0,lo,hi,job=None,slope=None,vtol=1e-9,maxev=40,max_bad=0):
    """Zero-current (floating) node: find V in [lo,hi] such that the total
    current into the device through the contacts `fidx` vanishes.
    The current is non-decreasing in V and, by the maximum principle, the root
    lies between the lowest and highest applied voltages (lo, hi), so the
    bracket needs no evaluation of its ends.  Safeguarded secant (Illinois
    regula falsi once both signs are known, bisection as a fallback); every
    evaluation continues from the previous solution.  max_bad > 0: give up
    after that many evaluations in a row fail to converge (past avalanche
    breakdown every one would).
    Returns (V, I_net, converged, evaluations, dI/dV estimate)."""
    a,b=float(lo),float(hi); fa=fb=None; side=0
    x=min(max(float(V0),a),b); hist=[]; ok=False; fx=0.; xr=x; found=False; bad=0
    for ev in range(1,maxev+1):
        if job is not None:
            if job.stop: break
            job.ev=ev
        fx,tol,ok,xr=_fl_eval(fidx,x)
        hist.append((xr,fx))
        if ok and abs(fx)<=tol: found=True; break
        bad=0 if ok else bad+1
        if max_bad and bad>=max_bad: break
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

def run_sweep(cons,tgt,Vs,Ve,Np,flt=None,job=None,on_point=None,vapp=None,refine=True,max_fail=0):
    """Sweep the node `tgt` (contact indices, driven together) from Vs to Ve in
    Np points, continuing from the present C state (setup_device first).
    Every point also records the electron and hole charge stored in the
    device (Qn, Qp in C: charge-control capacitances, the quasi-static f_T).
    vapp: applied voltage of every contact (default: the Con.V values); the
    non-swept, non-floating contacts are held there.
    flt: contact indices of a floating node solved for zero current at every
    point (e.g. an inverter output).  With refine=True an output jump larger
    than 8 % of the voltage span between neighbouring points is bisected
    (the steep part of a transfer curve gets its own points).
    on_point(k, rec) is called after every solved point (calling thread).
    max_fail: stop after this many points in a row fail to converge (at and
    just past avalanche breakdown the Gummel iteration fails at every point,
    and each failed point costs the solver its whole step-halving search);
    0 = never (the GUI uses 2 when impact ionisation is on).
    Returns a dict of numpy arrays ordered along the sweep (+ 'stop_fail':
    the voltage where max_fail stopped it, NaN if it did not)."""
    nc=len(cons); Np=max(1,int(Np))
    vapp=list(vapp) if vapp is not None else [c.V for c in cons]
    tgt=list(tgt); flt=[c for c in (flt or []) if c not in tgt]
    for c in range(nc):
        if c not in tgt and c not in flt: api_setV(c,vapp[c])
    Vv=np.linspace(Vs,Ve,Np) if Np>1 else np.array([float(Vs)])
    R=dict(V=[],I=[],F=[],conv=[],Vf=[],nev=[],t=[],Qn=[],Qp=[])
    fixed=[vapp[c] for c in range(nc) if c not in flt and c not in tgt]
    span=max(fixed+[Vs,Ve])-min(fixed+[Vs,Ve]) or 1.
    st={"slope":None}
    if job is not None: job.N=Np; job.extra=0
    def solve_at(Vt,guess):
        for c in tgt: api_setV(c,Vt)
        if flt:
            f=fixed+[Vt]; lo,hi=min(f),max(f)
            Vx,fx,ok,ne,sl=solve_floating(flt,min(max(guess,lo),hi),lo,hi,job,st["slope"],
                                          max_bad=1 if max_fail else 0)
            st["slope"]=sl
        else:
            ok=bool(api_resolve()); ne=1; Vx=float("nan")
        Vr=api_Vsol(tgt[0]); Qn,Qp=carrier_charge()
        rec=dict(V=(Vr if Vr==Vr else Vt),I=[api_conI(c) for c in range(nc)],
                 F=[api_conF(c) for c in range(nc)],conv=bool(ok),Vf=Vx,nev=ne,t=time.time(),Qn=Qn,Qp=Qp)
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
    prev=None; nfail=0; stop_fail=float("nan")
    for k,Vt in enumerate(Vv):
        if job is not None:
            if job.stop: break
            job.k=k; job.ev=0
        rec=solve_at(Vt,guess_at(Vt))
        nfail=0 if rec["conv"] else nfail+1
        if max_fail and nfail>=max_fail and k<len(Vv)-1:
            stop_fail=float(Vt); break
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
    out["stop_fail"]=stop_fail
    return out

# ── sweep analysis (GUI-independent; the I-V tab and validate.py use these) ──
def avalanche_bv(V,J,conv,thr=0.,vmin=3.,floor=None):
    """Breakdown voltage of a sweep solved with impact ionisation: where the
    current runs away.  V: swept voltage, J: |current density| (A/cm^2) of
    the observed terminal, conv: converged flags, floor: the current
    resolution (A/cm^2; currents below it are round-off and do not count) -
    all in sweep order.  Returns (BV, how) or (None, why).  Tried in turn:
      1. thr > 0: the current density reaches thr (the data-sheet way, BV
         at a stated leakage), interpolated in log J;
      2. the solver stops converging after the current has started to
         multiply (or the failed point already carries 100x the current):
         near breakdown (M ~ 10-50) the Gummel iteration fails.  BV is where
         1/J, extrapolated from the last two converged points, reaches zero
         (1/J ~ 1/M falls linearly to zero at BV - for a junction and for
         the open-base transistor alike), but at least the last converged
         voltage; without two resolved points, the middle of the interval;
      3. the current jumps tenfold between neighbouring converged points (a
         series resistance can keep the solver converging past breakdown).
    Only points with |V| >= vmin count (a forward-biased junction's current
    rises steeply too)."""
    V=np.asarray(V,float); J=np.abs(np.asarray(J,float)); cv=np.asarray(conv,bool)
    n=len(V)
    if n<3: return None,"too few points"
    res=np.isfinite(J)&(J>0)
    if floor is not None: res&=J>=np.asarray(floor,float)
    ok=cv&res; big=np.abs(V)>=vmin
    idx=np.where(ok)[0]
    if thr>0:
        for a,b in zip(idx[:-1],idx[1:]):
            if big[b] and J[a]<thr<=J[b]:
                la,lb=np.log(J[a]),np.log(J[b]); f=(np.log(thr)-la)/(lb-la) if lb>la else 1.
                return float(V[a]+f*(V[b]-V[a])),f"|J| reaches {thr:g} A/cm²"
    conv_idx=np.where(cv)[0]
    fail=[k for k in range(n) if not cv[k] and big[k] and np.any(conv_idx<k)]
    if fail:
        kf=fail[0]; prev=conv_idx[conv_idx<kf]
        k0=prev[-1]; kp=prev[-2] if len(prev)>=2 else None
        grow=bool(ok[k0]) and kp is not None and J[k0]>=1.3*(J[kp] if res[kp] else 0.)
        # or the failed point's (partial) state already carries a current far
        # above the last converged one (leakage below the resolution, then a
        # runaway - the open-base transistor)
        ref=J[k0] if res[k0] else (float(np.asarray(floor,float)[k0]) if floor is not None else 0.)
        jump=bool(np.isfinite(J[kf]) and J[kf]>0 and J[kf]>=100.*max(ref,1e-300))
        if grow or jump:
            BV=0.5*float(V[k0]+V[kf]); how=(f"the current runs away - the solver stops converging between "
                                            f"{V[k0]:.4g} and {V[kf]:.4g} V")
            if grow and res[kp] and J[k0]>J[kp] and V[k0]!=V[kp]:
                BV=float(V[k0])
                x=V[k0]+(V[k0]-V[kp])*(1./J[k0])/((1./J[kp])-(1./J[k0]))     # 1/J -> 0
                lim=V[k0]+3.*(V[k0]-V[kp])                  # at most three steps beyond the last point
                if (x-V[k0])*(V[k0]-V[kp])>=0 and abs(x-V[k0])<=abs(lim-V[k0]):
                    BV=float(x); how+=f"; 1/I extrapolates to zero at {x:.4g} V"
            return BV,how
    for a,b in zip(idx[:-1],idx[1:]):
        if big[b] and J[b]>=10.*J[a]:
            return float(V[a]),f"the current jumps {J[b]/J[a]:.0f}x between {V[a]:.4g} and {V[b]:.4g} V"
    return None,"no breakdown in the sweep range"

def early_voltage(V,I,conv):
    """Early voltage of an output characteristic (I_C vs V_CE at fixed V_BE,
    or I_D vs V_DS): a straight line I = a + b V fitted to the flat part of
    the curve - after the saturation knee, before any breakdown upturn -
    gives V_A = a/b (the line meets I = 0 at V = -V_A) and the output
    conductance g_o = b.  Works on magnitudes (PNP/PMOS curves too).
    Returns dict(VA, go, V0, V1, a, b, s) - s: the sign of V - or None."""
    V=np.asarray(V,float); I=np.asarray(I,float); cv=np.asarray(conv,bool)
    ok=cv&np.isfinite(I)&np.isfinite(V)
    if ok.sum()<5: return None
    V,I=V[ok],I[ok]
    sv=1. if V[np.argmax(np.abs(V))]>=0 else -1.
    m=(V*sv)>0
    if m.sum()<5: return None
    x=V[m]*sv; y=I[m]; o=np.argsort(x); x,y=x[o],y[o]
    if not np.all(np.diff(x)>0): x,iu=np.unique(x,return_index=True); y=y[iu]
    if len(x)<5: return None
    sy=1. if y[-1]>=0 else -1.; y=y*sy
    if y[-1]<=0: return None
    g=np.gradient(y,x)
    hi_=y>=0.3*y.max()                     # the flat part carries most of the current
    hi_[0]=hi_[-1]=False
    if not hi_.any(): return None
    j=int(np.argmin(np.where(hi_,g,np.inf)))
    gm=max(g[j],0.)
    tol=max(2.*gm,0.02*y[j]/max(x[j],1e-12))
    lo=hi=j
    while lo>0 and g[lo-1]<=tol: lo-=1
    while hi<len(g)-1 and g[hi+1]<=tol: hi+=1
    if hi-lo<2: return None
    b,a=np.polyfit(x[lo:hi+1],y[lo:hi+1],1)
    if not (b>0 and a>0): return None
    return dict(VA=float(a/b),go=float(b),V0=float(x[lo]*sv),V1=float(x[hi]*sv),a=float(a*sy),b=float(b*sy/sv),s=sv)

def bjt_gain_ft(V,Ib,Ic,Q,conv):
    """Current gain and quasi-static transit frequency along a Gummel sweep
    (V: the swept base voltage, V_CE held; Ib, Ic: base and collector
    currents; Q: the charge of the base's majority carrier stored in the
    device, C - holes for an NPN, electrons for a PNP).
       beta = I_C / I_B,   f_T = g_m / (2 pi C_in),
       g_m = dI_C/dV_BE,   C_in = dQ/dV_BE
    - the charge-control (low-frequency) estimate of the frequency where
    |h21| = 1: C_in holds the emitter and collector depletion charges and
    the stored minority charge, so the base pushing into the collector at
    high current (Kirk effect) shows as a fall of f_T past its peak.
    Returns dict(Ic, beta, fT, k) of numpy arrays over the converged points
    in |V| order (fT NaN where Q is unknown) - or None."""
    V=np.asarray(V,float); Ib=np.asarray(Ib,float); Ic=np.asarray(Ic,float)
    Q=np.asarray(Q,float) if Q is not None else np.full_like(V,np.nan)
    ok=np.asarray(conv,bool)&np.isfinite(Ic)&np.isfinite(Ib)&(np.abs(Ib)>0)
    if ok.sum()<3: return None
    V,Ib,Ic,Q=V[ok],Ib[ok],Ic[ok],Q[ok]
    o=np.argsort(np.abs(V)); V,Ib,Ic,Q=V[o],Ib[o],Ic[o],Q[o]
    keep=np.concatenate(([True],np.diff(V)!=0)); V,Ib,Ic,Q=V[keep],Ib[keep],Ic[keep],Q[keep]
    if len(V)<3: return None
    a=np.abs(Ic); beta=a/np.abs(Ib)
    # g_m = I_C d(ln I_C)/dV: no sinh-type error of a plain difference on coarse steps
    la=np.log(np.where(a>0,a,np.nan)); gm=np.abs(a*np.gradient(la,V))
    # f_T per interval between neighbouring points, dI_C/dQ = Delta I_C / Delta Q,
    # placed at the logarithmic mean current of the interval (IcT): on a coarse
    # Gummel sweep (I_C changing 2-3x per step) a point-wise value would sit at
    # the wrong current and shift the f_T peak (the Kirk current)
    if np.all(np.isfinite(Q)) and len(V)>=2:
        dI=np.abs(np.diff(a)); dQ=np.abs(np.diff(Q))
        fT=np.where(dQ>0,dI/(2.*np.pi*np.where(dQ>0,dQ,1.)),np.nan)
        a0,a1=a[:-1],a[1:]
        IcT=np.where((a0>0)&(a1>0)&(a0!=a1),(a1-a0)/np.log(np.where((a0>0)&(a1>0),a1/np.where(a0>0,a0,1.),2.)),0.5*(a0+a1))
    else: fT=np.full(max(len(V)-1,0),np.nan); IcT=0.5*(a[:-1]+a[1:])
    return dict(V=V,Ic=a,beta=beta,fT=fT,IcT=IcT,gm=gm)

KB_Q=8.617333262e-5          # Boltzmann constant / q (V/K)

def diode_ideality(V,J,T=300.,floor=None,decades=2.):
    """Ideality factor n and saturation current density J_0 of a forward
    diode branch J = J_0 exp(V/(n kT/q)): a straight line through ln J over
    the window of consecutive points (>= 3 points spanning >= `decades`
    decades, above 10x the resolution floor, V > 8 kT/q) with the smallest
    n - the ideal-diffusion part, away from depletion-region recombination
    at low bias and series resistance / high injection at high bias.
    Returns dict(n, J0, V0, V1) or None."""
    VT=KB_Q*(T if T else 300.)
    V=np.asarray(V,float); J=np.asarray(J,float)
    ok=np.isfinite(V)&np.isfinite(J)&(J>0)&(V>8.*VT)        # exp(V/nkT) >> 1 for n up to ~2
    if floor is not None: ok&=J>10.*np.abs(np.asarray(floor,float))
    V,J=V[ok],J[ok]
    if len(V)<3: return None
    o=np.argsort(V); V,J=V[o],J[o]; lj=np.log(J); need=decades*np.log(10.)
    best=None
    for i in range(len(V)-2):
        for k in range(i+2,len(V)):
            if lj[k]-lj[i]>=need: break
        else: break
        try: m,b=np.polyfit(V[i:k+1],lj[i:k+1],1)
        except Exception: continue
        if not m>0: continue
        n=1./(m*VT)
        if 0.5<n<10. and (best is None or n<best[0]): best=(n,float(np.exp(b)),float(V[i]),float(V[k]))
    if best is None: return None
    return dict(n=best[0],J0=best[1],V0=best[2],V1=best[3])

def subthreshold_swing(Vg,I,floor=None,decades=1.):
    """Subthreshold swing (mV/dec): the steepest straight line of log10|I|
    over windows of consecutive points spanning >= `decades` decades, all
    below 1 % of the largest current and above 10x the resolution floor (a
    single noisy pair near the floor cannot give a too-small SS).
    Returns dict(SS, V0, V1) or None."""
    Vg=np.asarray(Vg,float); I=np.abs(np.asarray(I,float))
    ok=np.isfinite(Vg)&np.isfinite(I)&(I>0)
    if floor is not None: ok&=I>10.*np.abs(np.asarray(floor,float))
    if ok.sum()<2: return None
    Imax=float(np.max(I[ok])); ok&=I<0.01*Imax
    Vg,I=Vg[ok],I[ok]
    if len(Vg)<2: return None
    o=np.argsort(Vg); Vg,I=Vg[o],I[o]; lg=np.log10(I)
    best=None
    for i in range(len(Vg)-1):
        for k in range(i+1,len(Vg)):
            if abs(lg[k]-lg[i])>=decades: break
        else: break
        if k-i>=2: m=np.polyfit(Vg[i:k+1],lg[i:k+1],1)[0]
        else: m=(lg[k]-lg[i])/(Vg[k]-Vg[i])
        if m==0: continue
        ss=abs(1e3/m)
        if best is None or ss<best[0]: best=(ss,float(Vg[i]),float(Vg[k]))
    if best is None: return None
    return dict(SS=best[0],V0=best[1],V1=best[2])

def threshold_voltage(Vg,I,Vds=None):
    """MOSFET threshold from a transfer curve (Ortiz-Conde et al. 2002):
    linear region (|V_DS| <= 0.1 V or unknown): tangent at maximum g_m,
    V_th = intercept - V_DS/2; saturation: tangent of sqrt|I_D| at its
    steepest point, V_th = intercept.  Returns dict(Vth, method, edge,
    gm, Vgm) or None; edge: the steepest point is at the end of the sweep."""
    Vg=np.asarray(Vg,float); I=np.abs(np.asarray(I,float))
    ok=np.isfinite(Vg)&np.isfinite(I)
    Vg,I=Vg[ok],I[ok]
    if len(Vg)<3: return None
    o=np.argsort(Vg); Vg,I=Vg[o],I[o]
    keep=np.concatenate(([True],np.diff(Vg)!=0)); Vg,I=Vg[keep],I[keep]
    if len(Vg)<3: return None
    gm=np.gradient(I,Vg)
    sat = Vds is not None and abs(Vds)>0.1
    if sat:
        y=np.sqrt(I); dy=np.gradient(y,Vg)
    else: y,dy=I,gm
    # the turn-on direction: the side where the current is larger
    sgn = 1. if I[-1]>=I[0] else -1.
    j=int(np.argmax(sgn*dy))
    if not sgn*dy[j]>0: return None
    Vint=Vg[j]-y[j]/dy[j]
    Vth = Vint if sat else Vint-(0.5*Vds if Vds is not None else 0.)
    return dict(Vth=float(Vth),method="sat" if sat else "lin",edge=j in (0,len(Vg)-1),
                gm=float(abs(gm[j])),Vgm=float(Vg[j]))

def contact_profile(s,cons,idx):
    """Current density through face contacts, node by node, from a pulled
    solution s: J into the device (A/cm^2) over each contact's nodes - how
    evenly a contact carries its current (emitter current crowding: the
    base current's lateral voltage drop debiases the middle of a wide
    emitter).  cons: the contacts of the solution, idx: the ones to show.
    Returns dict(parts=[(label, u, v, J(nu, nv), w(nu, nv))], axes=(U, V),
    skipped=[labels]) - u, v in nm along the face's two in-plane axes, w the
    node areas inside the contact (nm^2) - or None if no face contact."""
    x,y,z=s["x"],s["y"],s["z"]; Nx,Ny,Nz=len(x),len(y),len(z)
    comp={0:("Jx",1.),1:("Jx",-1.),2:("Jy",1.),3:("Jy",-1.),4:("Jz",1.),5:("Jz",-1.)}
    co={"X":x,"Y":y,"Z":z}
    def dual(c,lo,hi):                    # node widths clipped to [lo, hi]
        m=np.concatenate(([c[0]],0.5*(c[1:]+c[:-1]),[c[-1]]))
        return np.clip(np.minimum(m[1:],hi)-np.maximum(m[:-1],lo),0.,None)
    parts=[]; skipped=[]; axes=None
    for k in idx:
        c=cons[k]
        if c.face==6 or c.bc==2: skipped.append(c.label); continue
        i0,i1,j0,j1=c.to_idx(Nx,Ny,Nz,x,y,z)
        U,Vn=Con.FAXES[c.face]; cu,cv=co[U],co[Vn]
        key,sg=comp[c.face]; A=np.nan_to_num(s[key])*sg          # (Nz, Ny, Nx)
        e=(0 if c.face%2==0 else -1)
        if c.face in (0,1): J=A[j0:j1+1,i0:i1+1,e].T             # u = y, v = z
        elif c.face in (2,3): J=A[j0:j1+1,e,i0:i1+1].T           # u = x, v = z
        else: J=A[e,j0:j1+1,i0:i1+1].T                           # u = x, v = y
        u=cu[i0:i1+1]; v=cv[j0:j1+1]
        ulo,uhi=cu[0]+c.i0p*(cu[-1]-cu[0]),cu[0]+c.i1p*(cu[-1]-cu[0])
        vlo,vhi=cv[0]+c.j0p*(cv[-1]-cv[0]),cv[0]+c.j1p*(cv[-1]-cv[0])
        wu=dual(cu,ulo,uhi)[i0:i1+1]; wv=dual(cv,vlo,vhi)[j0:j1+1]
        if axes is None: axes=(U,Vn)
        elif axes!=(U,Vn): skipped.append(c.label); continue
        parts.append((c.label,u,v,J,np.outer(wu,wv)))
    if not parts: return None
    return dict(parts=parts,axes=axes,skipped=skipped)

def profile_stats(prof):
    """Mean (area-weighted), maximum, and centre / edge values of a contact
    profile; the edge and centre are read across the narrow direction of
    each contact (along u when the profile is one-dimensional)."""
    J=np.concatenate([p[3].ravel() for p in prof["parts"]]); W=np.concatenate([p[4].ravel() for p in prof["parts"]])
    W=np.where(W>0,W,0.)
    if W.sum()<=0: W=np.ones_like(J)
    mean=float((J*W).sum()/W.sum()); mx=float(J.max()); mn=float(J.min())
    ce=[]
    for _,u,v,Jp,Wp in prof["parts"]:
        alu=len(v)<3 or (len(u)>=3 and (u[-1]-u[0])<=(v[-1]-v[0]))     # across the narrow direction
        ax=1 if alu else 0
        w=Wp.sum(axis=ax); pr=(Jp*Wp).sum(axis=ax)/np.where(w>0,w,1.)
        if len(pr)>=3: ce.append((float(pr[len(pr)//2]),float(max(pr[0],pr[-1]))))
    return dict(mean=mean,max=mx,min=mn,ratio=mx/mean if mean else float("nan"),ce=ce)

api_getstate= _opt("api3_get_state",I,DP,DP,DP,default=0)
api_setstate= _opt("api3_set_state",I,DP,DP,DP,default=0)

def get_state():
    """(psi, phi_n, phi_p) of every node (V, flat arrays) of the C state, or None."""
    N=_PRISM.N if _PRISM is not None else api_Nx()*api_Ny()*api_Nz()
    bufs=[(D*N)() for _ in range(3)]
    if not api_getstate(*bufs): return None
    return [np.frombuffer(b,dtype=np.float64,count=N).copy() for b in bufs]

def set_state(state):
    """Make (psi, phi_n, phi_p) the solution at the present contact voltages."""
    arrs=[np.ascontiguousarray(a,dtype=np.float64) for a in state]
    return bool(api_setstate(*[a.ctypes.data_as(DP) for a in arrs]))

def _interp_axis(A,c_old,c_new,axis):
    c_old=np.asarray(c_old,float); c_new=np.asarray(c_new,float)
    i=np.clip(np.searchsorted(c_old,c_new)-1,0,len(c_old)-2)
    w=np.clip((c_new-c_old[i])/(c_old[i+1]-c_old[i]),0.,1.)
    shp=[1]*A.ndim; shp[axis]=len(c_new); w=w.reshape(shp)
    return np.take(A,i,axis=axis)*(1.-w)+np.take(A,i+1,axis=axis)*w

def transfer_state(state,old,new):
    """A tensor-mesh state on the lines old=(xs,ys,zs) interpolated (linearly)
    onto the lines new."""
    (xo,yo,zo),(xn,yn,zn)=old,new; out=[]
    for a in state:
        A=np.asarray(a).reshape(len(zo),len(yo),len(xo))
        A=_interp_axis(A,xo,xn,2); A=_interp_axis(A,yo,yn,1); A=_interp_axis(A,zo,zn,0)
        out.append(np.ascontiguousarray(A.ravel()))
    return out

AMR_HMIN=0.5        # nm: no interval is split below twice this

def amr_marks(lines,state,s,kappa=1.,T=300.,mat=None):
    """Where a tensor-mesh solution is under-resolved.  Indicator at each node
    and axis: how far a potential deviates from the straight line through
    the node's two neighbours (curvature x h^2), in units of kT/q - the
    electrostatic potential, the band edge each carrier moves in (Ec, Ev:
    with band-gap narrowing and the quantum-confinement potential) and the
    quasi-Fermi potentials, the carrier quantities only where that carrier
    is within 2 decades of the majority; semiconductor nodes of one material
    only (mat: material index per node).  Unresolved space charge (h above
    ~ the Debye length), inversion layers and current crowding show up; a
    uniform field (a depleted drift region, a resistor) does not.  Returns
    [(eta, axis, [new positions nm], cost in nodes)] for the intervals next
    to nodes with eta > kappa, largest eta first; an interval is split into
    sqrt(eta/kappa) parts (the deviation falls as h^2), at most 8."""
    VT=8.617333262e-5*T; nz,ny,nx=len(lines[2]),len(lines[1]),len(lines[0])
    shape=(nz,ny,nx)
    psi,phn,php=(np.asarray(a).reshape(shape) for a in state)
    semi=np.asarray(s["semi"],bool).reshape(shape)
    n=np.nan_to_num(np.asarray(s["n"],float)).reshape(shape); p=np.nan_to_num(np.asarray(s["p"],float)).reshape(shape)
    Ec=np.nan_to_num(np.asarray(s["Ec"],float)).reshape(shape); Ev=np.nan_to_num(np.asarray(s["Ev"],float)).reshape(shape)
    big=np.maximum(n,p); sn=semi&(n>=1e-2*big); sp=semi&(p>=1e-2*big)
    mid=np.zeros(shape,np.int32) if mat is None else np.asarray(mat).reshape(shape)
    out=[]; cnt=(nx,ny,nz)
    for q,axn in enumerate("xyz"):
        c=np.asarray(lines[q],float); m=len(c)
        if m<3: continue
        ax=2-q; h=np.diff(c); hm=h[:-1]; hp=h[1:]
        bs=[1,1,1]; bs[ax]=m-2; hm=hm.reshape(bs); hp=hp.reshape(bs)
        tk=lambda F,a,b: np.take(F,np.arange(a,b),axis=ax)
        def dev(F):
            return np.abs(tk(F,1,m-1)-(hp*tk(F,0,m-2)+hm*tk(F,2,m))/(hm+hp))
        same=(tk(mid,0,m-2)==tk(mid,1,m-1))&(tk(mid,1,m-1)==tk(mid,2,m))
        def all3(M): return tk(M,0,m-2)&tk(M,1,m-1)&tk(M,2,m)&same
        e=np.where(all3(semi),dev(psi),0.)
        mn_=all3(sn); mp_=all3(sp)
        e=np.maximum(e,np.where(mn_,np.maximum(dev(phn),dev(Ec)),0.))
        e=np.maximum(e,np.where(mp_,np.maximum(dev(php),dev(Ev)),0.))
        oth=tuple(a for a in range(3) if a!=ax)
        en=e.max(axis=oth)/VT                       # interior nodes 1..m-2
        ei=np.zeros(m-1); ei[:-1]=np.maximum(ei[:-1],en); ei[1:]=np.maximum(ei[1:],en)
        cost=int(np.prod([cnt[k] for k in range(3) if k!=q]))
        for j in np.nonzero(ei>kappa)[0]:
            hn=h[j]*1e9
            if hn<2.*AMR_HMIN: continue
            ns=int(min(8,np.ceil(np.sqrt(ei[j]/kappa)),max(2,hn//AMR_HMIN)))
            pos=[(c[j]+(c[j+1]-c[j])*f/ns)*1e9 for f in range(1,ns)]
            out.append((float(ei[j]),axn,pos,cost*len(pos)))
    out.sort(key=lambda r:-r[0])
    return out

def amr_run(t,cons,flt=None,job=None,vapp=None,solver=(2000,1000,500,1e-8,1e-6,1.,1.),models=7,T=300.,
            graded=True,passes=3,max_nodes=None,kappa=1.,on_pass=None):
    """Solve at the applied voltages (the device set up from t), then refine
    the rectangular mesh where the solution is under-resolved (amr_marks) and
    solve again from the interpolated solution, up to `passes` times or until
    nothing is under-resolved or the node budget max_nodes (default: twice
    the starting mesh) is used.  t['amr'] collects the added lines (nm).
    Returns (ok, V_float, [(lines added, nodes, converged) per pass])."""
    vapp=list(vapp) if vapp is not None else [c.V for c in cons]
    ok,Vx=run_point(cons,flt,job,vapp)
    info=[]
    if t.get("mesh")=="tri" or _PRISM is not None: return ok,Vx,info
    lines=tuple(np.asarray(a,float) for a in make_grid(t,graded))
    N0=len(lines[0])*len(lines[1])*len(lines[2])
    cap=int(max_nodes) if max_nodes else 2*N0
    for ps in range(int(passes)):
        if (job is not None and job.stop) or not ok: break
        s=pull(); st=get_state()
        if st is None: break
        own=owner_grid(t["regs"],*lines)
        mats=np.array([MATS[r.mat] for r in t["regs"]]+[-1])[np.where(own<0,len(t["regs"]),own)]
        cand=amr_marks(lines,st,s,kappa,T,mats); sel=[]
        cnt=[len(lines[0]),len(lines[1]),len(lines[2])]            # lines per axis as they grow
        for eta,axn,pos,_ in cand:
            q="xyz".index(axn); trial=list(cnt); trial[q]+=len(pos)
            if trial[0]*trial[1]*trial[2]>cap: continue
            sel.append((axn,pos)); cnt=trial
        if not sel: break
        amr=t.get("amr") or {}; t["amr"]=amr
        for axn,pos in sel: amr.setdefault(axn,[]).extend(pos)
        if flt and Vx is not None:
            for c in flt: vapp[c]=Vx
        volt={cons[c].label:vapp[c] for c in range(len(cons))}
        new=setup_device(t,solver=solver,models=models,T=T,graded=graded,volt=volt)
        new=tuple(np.asarray(a,float) for a in new)
        set_state(transfer_state(st,lines,new)); lines=new
        ok,Vx=run_point(cons,flt,job,vapp)
        if not ok and not (job is not None and job.stop):
            # the interpolated start did not converge: the refined device from equilibrium
            setup_device(t,solver=solver,models=models,T=T,graded=graded,volt=volt)
            ok,Vx=run_point(cons,flt,job,vapp)
        info.append((sum(len(p) for _,p in sel),len(lines[0])*len(lines[1])*len(lines[2]),bool(ok)))
        if on_pass: on_pass(ps,info[-1])
    return ok,Vx,info

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
    cons=[(c.face,c.i0p,c.i1p,c.j0p,c.j1p,c.k0p,c.k1p,c.bc,c.pm,c.tox,c.nox,repr(c.stack),_shape_key(c),
           getattr(c,"metal","")) for c in t["cons"]]
    shs=[(s.axis,s.pos,s.a0,s.a1,s.b0,s.b1,s.sigma) for s in t.get("sheets",[])]
    return repr((t["Lx"],t["Ly"],t["Lz"],t["Nx"],t["Ny"],t["Nz"],regs,cons,shs,
                 sorted(t.get("dense_x",[])),sorted(t.get("dense_y",[])),sorted(t.get("dense_z",[])),
                 repr(t.get("nodes")),repr(t.get("clear")),repr(t.get("refine")),tuple(solver),models,T,graded,
                 t.get("mesh","rect"),t.get("prism"),t.get("extrude"),
                 [d.key() for d in (t.get("defects") or [])],t.get("tau0"),
                 repr(sorted((k,sorted(v)) for k,v in (t.get("amr") or {}).items()))))

# ══════════════════════════════════════════════════════════════
#  GUI TOOLKIT: fonts, two-way scrolling pages, tooltips, help text
# ══════════════════════════════════════════════════════════════
import queue, json, copy, textwrap
import tkinter.font as tkfont
from matplotlib.ticker import FuncFormatter

VERSION="0.2"
SETTINGS_FILE=os.path.join(os.path.expanduser("~"),".semisim3d.json")
# label-based buttons (FlatButton) and themed check boxes: macOS draws the
# classic Tk ones natively, without the app's colours; SEMISIM_FLAT_BUTTONS=1
# uses them elsewhere too (tests)
FLAT_BUTTONS=IS_MAC or os.environ.get("SEMISIM_FLAT_BUTTONS")=="1"
MOD="Cmd" if IS_MAC else "Ctrl"           # shortcut modifier shown in menus and help
def _aqua(w):
    """Is Tk drawing with macOS's native toolkit?  (SEMISIM_FAKE_AQUA=1 takes
    the macOS code paths elsewhere, for tests; the Mac-only Tk commands then
    fail quietly.)"""
    if os.environ.get("SEMISIM_FAKE_AQUA")=="1": return True
    try: return w.tk.call("tk","windowingsystem")=="aqua"
    except Exception: return False
def _x11(w):
    try: return w.tk.call("tk","windowingsystem")=="x11"
    except Exception: return False
ARROWS_2D=(4,120)          # Current tab: arrows across the in-plane map (min, max)
ARROWS_3D=(2,40)           # ... along x in the 3-D field (y: 3/4, z: 5/8 of it)
ARROWS_3D_MAX=4000         # total 3-D arrows (mplot3d redraws every arrow on each rotation step)

def _lin_axis(F,c,g,axis):
    """Linear interpolation of the array F along `axis` (node coordinates c,
    increasing) at the points g.  Where a neighbour is NaN (no carriers: an
    oxide, a void) the nearer node's value is used."""
    c=np.asarray(c,float); g=np.asarray(g,float); n=len(c)
    if n<2: return np.take(F,np.zeros(len(g),int),axis=axis)
    i=np.clip(np.searchsorted(c,g)-1,0,n-2)
    d=c[i+1]-c[i]; t=np.clip((g-c[i])/np.where(d>0,d,1.),0.,1.)
    A=np.take(F,i,axis=axis); B=np.take(F,i+1,axis=axis)
    sh=[1]*F.ndim; sh[axis]=len(g); t=t.reshape(sh)
    with np.errstate(invalid="ignore"):
        out=(1.-t)*A+t*B
    bad=~np.isfinite(out)
    if bad.any(): out=np.where(bad,np.where(t>0.5,B,A),out)
    return out

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
            if IS_MAC:          # a help-style window that never takes the focus
                try: self.tw.tk.call("::tk::unsupported::MacWindowStyle","style",self.tw._w,"help","noActivates")
                except tk.TclError: pass
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

class FlatButton(tk.Label):
    """Push button that keeps its colours on every platform.  macOS draws
    tk.Button natively and ignores -background, so the white captions of
    the coloured buttons would sit on a white button.  A Label has the same
    -state, -activebackground and -disabledforeground options; this adds the
    -command, the hover highlight, invoke() and the space/enter keys."""
    def __init__(self,master=None,command=None,**kw):
        kw.pop("relief",None); kw.pop("overrelief",None)
        kw.setdefault("cursor","hand2"); kw.setdefault("takefocus",1)
        self._cmd=command
        super().__init__(master,relief="flat",**kw)
        self._hover=False
        self.bind("<Enter>",self._enter,add="+"); self.bind("<Leave>",self._leave,add="+")
        self.bind("<ButtonRelease-1>",self._release,add="+")
        for k in ("<space>","<Return>","<KP_Enter>"): self.bind(k,lambda e:self.invoke(),add="+")
    def _enter(self,_=None):
        self._hover=True
        if str(super().cget("state"))=="normal": super().configure(state="active")
    def _leave(self,_=None):
        self._hover=False
        if str(super().cget("state"))=="active": super().configure(state="normal")
    def _release(self,e):
        if 0<=e.x<self.winfo_width() and 0<=e.y<self.winfo_height(): self.invoke()
    def invoke(self):
        if str(super().cget("state"))!="disabled" and self._cmd: return self._cmd()
    def configure(self,cnf=None,**kw):
        if isinstance(cnf,dict): kw={**cnf,**kw}; cnf=None
        elif isinstance(cnf,str): return self._cmd if cnf=="command" else super().configure(cnf)
        if "command" in kw: self._cmd=kw.pop("command")
        kw.pop("relief",None); kw.pop("overrelief",None)
        if kw.get("state")=="normal" and self._hover: kw["state"]="active"
        return super().configure(**kw)
    config=configure
    def cget(self,key):
        if key=="command": return self._cmd
        v=super().cget(key)
        return "normal" if key=="state" and str(v)=="active" else v
    __getitem__=cget
    def __setitem__(self,key,val): self.configure(**{key:val})

def render_markup(txt,body):
    """Tiny markup -> Tk text tags:  '# ' heading, '## ' / '### ' subheading,
    '- ' bullet, '  - ' sub-bullet, '  ' (two spaces) preformatted line,
    `code` spans."""
    import re
    for line in body.strip("\n").split("\n"):
        if line.startswith("### "):
            txt.insert("end",line[4:]+"\n","h2"); continue
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
("Quick start", f"""
# Quick start
- Pick a device from the `Device ▾` menu in the toolbar (or `Devices` in the menu bar). The structure appears in the `Structure` tab.
- Press `▶ Run` (F5) to solve the device at the voltages listed in the editor's `Contacts` page.
- Press `↔ Sweep` (F6) to run the I-V sweep defined in the `Sweep` page. The curve builds up live in the `I-V` tab.
- Press `■ Stop` (Esc) at any time. The solver keeps the last converged bias point, so all plots stay valid.
- Open or close the editor panel with `☰` ({MOD}+E). With the panel closed the plots use the whole window.
- F11 makes the window fill the screen (and back). On WSL use F11 rather than the window's maximise button - WSLg mis-routes clicks in maximised windows (the maximise button is turned into F11 automatically).
- Every panel scrolls both ways. Use the scroll bars, the mouse wheel, or Shift + wheel for sideways scrolling.
- Text too small or too large? Use View → Text size, or {MOD} + / {MOD} − ({MOD} 0 resets). Plot labels scale too.

Typical workflow for a new design: load the closest template, change regions and contacts, press Run, inspect the Bands/Current/Field tabs, then sweep.
"""),
("Status and terms", """
# Status and terms
SemiSim 3D 0.2 is an early public version - a test version: it is still being tested, so read its results with care.
- The physics is checked against textbook theory (`validate.py`, 213 checks), not yet against measurements of real devices.
- Results can be wrong, and some effects are not modelled at all (see Limitations). Check anything important against measurements, the literature or an established simulator.
- The `RW:` templates are modelled on the device class and ratings of real parts, not on their real dies; currents are per simulated cell.
- macOS support and the WSL2 full-screen workaround are new and have been tried on few computers so far.
- The format of saved device files may change between versions.
- Problems and suggestions: https://github.com/kotmariusz/Semiconductor_devices_simulator/issues - say what you did, attach the saved device (File → Save device) and the output of `python3 semibuild.py --check`.
## Licence
SemiSim 3D is released under the MIT License (the LICENSE file in the SemiSim 3D folder): free to use, copy, change and share, for any purpose, as long as the copyright and licence notice stay with it. There is no warranty of any kind: whoever uses it is responsible for how they use its results.
Part numbers and company names in the `RW:` templates only identify the parts the templates are modelled on; they belong to their owners, who are not connected with this project.
## Sources
The models, the material and defect parameters, the device structures and the test references come from the literature listed under Sources (also `docs/SOURCES.md`).
"""),
("Window layout", """
# Window layout
## Toolbar
- `☰` shows or hides the editor panel (drawer). Drag the divider to resize it.
- `Device ▾` opens the template library by category. `Browse…` shows full descriptions.
- `▶ Run`, `↔ Sweep` and `■ Stop` control the solver. The progress bar shows the bias ramp (Run) or the sweep points (Sweep). The text beside it shows the continuation step, the Gummel iteration and the residual.
## Editor panel (drawer) pages
- `Device`: the template description, domain size, base mesh, mesh type (rectangular or triangular prisms), adaptive mesh refinement and mesh statistics.
- `Regions`: the material/doping blocks. Later rows override earlier ones where they overlap.
- `Contacts`: electrodes, their voltages, work functions, gate oxides and nets.
- `Defects`: voids, particles, contamination, radiation damage, interface traps, dislocations - placed by hand or scattered at random (see Defects).
- `Solver`: temperature, physical models, background SRH lifetime, iteration limits, tolerances and CPU threads.
- `Sweep`: the swept terminal, the voltage range, a floating (zero-current) terminal and the plotted current.
## Plot tabs
- `Structure`: 3-D view and a cross-section with the real mesh (lines, or the triangles of a prism mesh). Pick the plane XY/YZ/XZ and its position.
- `3D Slicer`: any solution quantity on three orthogonal planes.
- `Bands`: Ec, Ev and the quasi-Fermi levels along a line cut, plus n and p.
- `Current`: current density along a cut, recombination, the in-plane current field and a 3-D arrow field (arrow counts adjustable).
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
The I-V tab reports the subthreshold swing SS, V_th (linear extrapolation at maximum g_m minus V_DS/2, or √I_D extrapolation in saturation) and I_on/I_off; the semi-log panel also shows the gate tunnelling current.
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
## Adaptive mesh refinement (Device page)
With `Refine the mesh to the solution on ▶ Run` ticked (rectangular mesh), every Run is followed by a check of the solution: at each node and along each axis, how far the electrostatic potential, the band edges the carriers move in (with band-gap narrowing and the quantum-confinement potential) and the quasi-Fermi potentials of the carriers that matter there bend away from the straight line through the two neighbours - a measure of how much detail falls inside one cell. Where that exceeds the threshold (in kT/q, 26 mV at 300 K), the mesh does not resolve the solution: an unresolved depletion edge, an inversion layer, current crowding at a junction. Those intervals are split - into √(deviation/threshold) parts, since the deviation falls with the square of the spacing - and the device is solved again, starting from the old solution interpolated onto the new mesh (if that does not converge, from equilibrium).
- `Passes`: how many check-and-refine rounds (default 3).
- `Max nodes ×`: the budget, as a multiple of the unrefined node count (default 2); the largest deviations are refined first.
- `Threshold`: lower = finer (default 1 kT/q).
- A uniform field (a depleted drift region), a resistor or a resolved junction needs no refinement and gets none.
- On a rectangular mesh a new line runs through the whole device, so the budget counts whole planes.
- The added lines stay for sweeps and later runs until `Clear refinement` or another template; they are saved with the device. Refine at the bias that matters most (the highest current, the region you study), then sweep.
- Quantities that hinge on something the mesh cannot fix are not improved: a model discontinuity (the quantum-confinement potential ends abruptly where a face gate ends), or a part of the device the indicator does not flag. When in doubt, compare with a run on a uniformly finer mesh.
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
## Contact metal
`Metal` makes the contact out of a real metal or silicide: Al, Ti, Ag, Cr, W, Mo, Cu, Au, Pd, Ni, Pt, TiN, TaN, TiSi2, CoSi2, NiSi, PtSi, ErSi2, n+ or p+ poly-Si. The line under the field shows the Schottky barrier the metal forms on the semiconductor under the contact, and where the number comes from:
- measured: barrier heights on n-type Si, GaAs, 4H-SiC and GaN from the literature (Sze & Ng; the SiC/GaN contact literature), for example Al/n-Si 0.72 eV, Ti/n-Si 0.50, PtSi/n-Si 0.84, Ni/4H-SiC 1.5 eV;
- pinning model: any other pair, φBn = S(φm − χ) + (1 − S)φpin (Cowley-Sze; S = 0.27 and φpin = 0.80 eV on Si; S ≈ 0.07 on GaAs - strongly pinned; S ≈ 0.4 on SiC and GaN).
The hole barrier is Eg − φBn. On a semiconductor face, both `Ohmic` and `Schottky` contacts with a metal become physical metal contacts (Solver page, metal-contact physics): the metal's barrier sits at the interface, and electrons and holes cross it by thermionic emission over the image-force-lowered barrier and by tunnelling through it. The doping then decides what the contact does: on 1e20 cm⁻³ Si the barrier is ~3 nm thin and electrons tunnel straight through - an ohmic contact (ρc ≈ 1e-7..1e-6 Ω·cm²); on 1e16 it is a Schottky diode. A `Gate` with a metal takes the metal's work function as φm. A Box contact keeps an ideal ohmic contact (Ohmic) or the fixed barrier (Schottky). `(ideal / custom φm)` keeps the classic contacts: ideal ohmic, or a Schottky barrier from the φm you type. Real barriers vary by ±0.1 eV with surface cleaning and anneal: type φm to override.
- `Gate stack` (face gates): the layers between the semiconductor and the metal, semiconductor side first, in nm, e.g. `SiO2 0.5, HfO2 2.5`. The gate capacitance uses the stack's EOT (tox is then computed) and gate tunnelling sees every layer. Empty means one SiO2 layer of thickness tox. A Box gate needs no stack: tunnelling follows the meshed insulator regions between the semiconductor and the metal.
Common work functions: n+ poly 4.05–4.1 eV, p+ poly 5.1–5.2 eV, TiN 4.5–4.6 eV, Ni 5.1 eV.
## Nets
Contacts with the same `Net` name are wired together outside the simulated cell. They always carry the same voltage, they are swept together, and their currents add up in the plots. Editing the voltage of one member updates the whole net.
Examples: the NMOS and PMOS gates of an inverter (`Vin`), the three sides of a FinFET gate (`Gate`), and a source and substrate contact tied to ground.
## Sign convention
A terminal current is positive when current flows INTO the device through that terminal. The currents of all terminals add up to zero (Kirchhoff).
"""),
("Defects", """
# Defects
The `Defects` page puts defects into any device. A defect is a physical effect from the library (literature values, every number editable) placed in a shape.
## What a defect does
- Deep levels (traps): Shockley-Read-Hall recombination through the level - any energy, its own capture cross sections σn and σp, thermal velocities from the effective masses - and the charge the level holds. An acceptor is negative when it holds an electron, a donor positive when it is empty; the occupancy follows the local electron and hole densities (steady state) and enters Poisson's equation, so traps compensate doping, deplete regions and shift thresholds. A neutral level only recombines.
- Interface traps, D_it (cm⁻²eV⁻¹): a constant density across the gap at the semiconductor/insulator interfaces inside the shape and under face gates, as 24 levels - acceptor-like above midgap, donor-like below. They add C_it = q²D_it to the subthreshold swing and shift the threshold by the charge they hold.
- Fixed charge: at interfaces (q/cm², e.g. oxide charge from processing or radiation) or in a volume (q/cm³).
- Inclusions: the shape becomes another material - a void (ε = 1), SiO2 (an oxide precipitate) or a metal particle (ε → ∞, a floating equipotential) - with interface states on its surface. A metal particle's states (5e13 cm⁻² each side of the pinning level) pin the Fermi level at the metal's Schottky barrier on that semiconductor: it depletes its surroundings like a floating Schottky contact and recombines like one.
- Extended defects: dislocations (deep levels along a line, states per cm), stacking faults and grain boundaries (states per cm² of the plane, or a D_it band).
- Radiation damage: neutron fluence (three-level Perugia model, level densities = η × Φ) and total ionizing dose (oxide charge + D_it).
## Library
- Voids & inclusions: void, oxide precipitate (SiO2), metal/silicide particle (any metal of the Contacts page), Cu3Si precipitate (copper contamination), COP (crystal-originated void), micropipe (4H-SiC).
- Point defects & contamination: a custom deep trap, interstitial iron Fe_i, the FeB pair, gold and platinum (lifetime control).
- Radiation damage: neutron fluence (Perugia), ionizing dose.
- Interfaces: Si/SiO2 D_it, 4H-SiC/SiO2 (D_it + near-interface acceptor traps), oxide fixed charge, volume charge.
- Extended defects: dislocation (Si), threading dislocation (GaN), stacking fault (4H-SiC), grain boundary (poly-Si).
- SiC · GaN · GaAs: Z1/2 and EH6/7 (SiC lifetime killers), C_N and Fe_Ga (semi-insulating GaN buffers), EL2 (GaAs).
The description under the selection gives each value and its source. Measured capture cross sections scatter by up to 10x (much more for Pt), so treat them as starting points.
## Shapes and placement
- `sphere` (diameter), `box` (dx, dy, dz), `cylinder` (diameter, and length along its axis - 0 = through the device), `slab` (thickness along its normal axis and two widths - 0 = across the device), `everywhere` (the whole host material).
- Centre in nm; `centre` puts it in the middle of the device. `Host`: the semiconductor the defect lives in (`any` = every semiconductor). `Level`: type, energy (below Ec, above Ev or from midgap Ei) and σn/σp; editing it replaces the first level of the preset.
- Manual: fill the form, `+ Add`; select a row to load it, `✎ Update`. Later rows override earlier ones where inclusions overlap.
- Random: `Random scatter` places copies of the form's defect at random centres in a box - a number, or a density (per cm³ of the box; per cm² for lines through the device, e.g. a threading-dislocation density), a seed (the same seed gives the same layout) and a size spread. The copies form a group that `Remove group` deletes at once.
## Mesh and exactness
The number of states is exact on any mesh: deep levels in a shape smaller than the cells are spread over the nearest nodes, a dislocation over the nearest row of nodes, and interface states carry the true surface area of the shape (the cell faces around a meshed sphere add up to ~1.5x its area). An inclusion replaces whole control volumes, so the graded mesh is made denser through the first 40 inclusions, lines and planes (boxes get exact faces). An inclusion smaller than the cells around it would grow several times its size, so it acts through its surface states only; the Defects page says so after a run, together with inclusions meshed with a very different volume. Raise Nx/Ny/Nz (or add dense planes) to resolve them.
## Background lifetime
`Background SRH τ0` (Solver page) is the lifetime of the semiconductors' built-in midgap recombination: empty = material defaults (Si 10 µs), `1e-3` = float-zone-quality silicon, `none` = no background SRH at all - only the defects recombine. With the 10 µs default a weak defect can drown in the background.
## Seeing them
- `Structure`: defects cut by the section are filled, the others dotted (their projection on the plane); the 3-D view shows all of them; hovering names the defect under the cursor.
- Solution maps (3D Slicer, Current) outline the defects cut by the plane when `Regions` is on. Slicer quantities `trap density` (all deep levels, cm⁻³) and `trapped charge`.
- Defects are saved with the device (File → Save device).
## Templates
- `Defects: leaky p-n junction (particle, void, dislocation)`: a NiSi precipitate in the depletion region raises the reverse leakage and the low-bias recombination current by three orders of magnitude; the void and the dislocation add their own recombination.
- `Defects: NMOS after irradiation (oxide charge + D_it)`: the threshold shift and swing of a total-ionizing-dose damaged MOSFET.
## Validation (validate.py, section 15)
Compensation of n-Si by a deep acceptor against charge neutrality (exact); the forward current of a p+/n diode whose lifetime comes from a trap level against Shockley's long-base theory (within 0.7 %, also with no background SRH); the swing increase from D_it against kT/q ln10 C_it/C_ox (1 %) and the threshold shift from oxide charge against −qQ/C_ox (0.02 %) on the NMOS; a void defect against the same box as a region; a NiSi particle pinning the Fermi level at 0.65 eV; and the exact number of states of lines, spheres and void surfaces on both mesh types, resolved and below the mesh size.
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
- `BV knee`: a current density (A/cm²) for breakdown. With impact ionisation on (Solver page) the breakdown voltage is where the current reaches it - 0: where the current runs away. With impact ionisation off it only places a heuristic knee marker.
- With impact ionisation on, a sweep ends once two points in a row fail to converge (at and just past avalanche breakdown every point would fail, and each failed point costs the solver its whole step-halving search).
## Stop (Esc)
Stops within one Gummel iteration. A sweep keeps the points already solved, and the device keeps its last converged state.
## Extracted parameters
The I-V tab and the Results tab report what fits the sweep:
- diode ideality factor and J0 (forward exponential branch);
- subthreshold swing, V_th and I_on/I_off (gate sweeps);
- V_M, gain and noise margins (transfer curves);
- the Early voltage V_A = I/g_o − V of output curves (collector or drain swept and observed): a straight line through the flat part, drawn dotted;
- for a Gummel plot (base swept, collector observed): the current gain β = I_C/I_B and the quasi-static f_T against I_C (the BJT panel), f_T,max and where β halves;
- the avalanche breakdown voltage when impact ionisation is on (else a heuristic current knee).
Pinned sweeps get the same parameters, listed below those of the last sweep; `Table…` puts them side by side with their differences.
## Comparing designs
- `Pin sweep` (I-V tab) keeps the last sweep as a reference: later sweeps of the same swept terminal are drawn over it, dashed in the pin's colour, on the linear, semi-log and BJT panels.
- `Pins…`: show or hide each pin (click 'Shown'), rename it (double-click), change its colour, show the fields of its last point in the solution tabs ('Show fields'), export all pins to CSV, remove them.
- `A/B compare…`: capture two designs and run the Sweep page's sweep on both. A design is the whole editor state - structure, mesh, contact voltages, physics models, temperature and an open (floating) terminal - captured with 'Capture the editor's device'; 'Load into the editor' brings it back. Templates of the 'BJT effects & breakdown' group have variant pairs: pick one under 'Variants', press '→ A', pick the other, '→ B', then '▶ Run A and B'. The two results become pins (A blue, B pink); the older pins are hidden.
- A design B can differ from A in anything: a doping, a finger layout, one physics switch (impact ionisation on/off), the temperature.
- `Table…` (I-V tab): the parameters of the shown pins and the last sweep side by side - device, sweep, temperature, physics models ('default' = the Solver page's default switches, + and − what differs), then V_th (max-g_m and constant current), SS, g_m, I_on, I_off, diode n and J_0, V_A, BV, β,max, f_T,max, V_M, V_OH/V_OL, V_IL/V_IH, noise margins and the output current at the sweep end. Choose the `Reference` column: the others show their difference from it - voltages in mV or V, currents and frequencies as ratios (×), gain, f_T, V_A and BV also in per cent. `Show`: values, values and differences, or differences only; `Hidden pins too` adds hidden pins. `Copy` (tab-separated, for a spreadsheet) and `Export CSV…` (numbers in V, A, S, Hz). `*` = a value with a note in the I-V tab, e.g. still rising at the sweep end.
"""),
("Plots", """
# Plots
- `Structure`: choose the cross-section plane (XY, YZ, XZ) and slide its position. Mesh lines are the real nodes. Hover over the section to read the region, material and doping under the cursor (and the defect, if any). Defects cut by the plane are filled, the others dotted.
- Every cut control has an axis selector and two position sliders for the other two axes. The label shows the position in nm.
- `Bands`: Ec and Ev (solid), the quasi-Fermi levels Efn and Efp (dashed) and a combined Ef, weighted towards the majority carrier. The energy zero is the vacuum level at φ = 0, so absolute values depend on the gauge but differences are physical. The n and p overlay uses a log scale on the right axis.
- `3D Slicer`: pick a quantity and slide the plane. The three small maps are the XY, YZ and XZ planes through the slice point, drawn on the real (non-uniform) mesh.
- `Current`: the in-plane map shows |J| with direction arrows. Arrow length is log-scaled so weak and strong regions are both visible. The number boxes next to `2-D arrows` and `3-D arrows` set how many arrows are drawn: 2-D, the arrows across the map (4-120; the rows follow the map's shape); 3-D, the arrows along x (2-40; depth y gets 3/4 and z 5/8 of them, never more than the mesh lines along an axis, at most 4000 arrows). The current is interpolated between the mesh nodes, so the count is not limited by the mesh; denser arrows are drawn shorter and thinner.
- `Panels ▾` (every plot tab): hide panels to give the others the room - the Structure tab's 3-D view, cross-section and legend; the Slicer's 3-D view and its three planes; the band diagram and ribbon; the five Current panels; E-field and mobility; the linear, semi-log and BJT panels of the I-V tab.
- `Layers…` (Structure and Current tabs): click 'Shown' to hide or show a region, contact or the defects in the 3-D view. Select one or more parts (Ctrl/Shift for several) and move 'Opacity of the selected parts' to give them their own opacity - a see-through substrate under an opaque channel, say; 'Auto' hands them back to the sliders for all regions and all contacts ('See-through contacts' sets 0.15). 'Insulators in 3-D' draws the oxides too; 'Apply to the cross-section too' uses the same choices in the section.
- `Current` → `Map`: what lies under the arrows of the in-plane map - `|J|`, the device `structure` (region colours) or `|J| + structure`. `Structure in 3-D` draws the regions and contacts around the 3-D arrows. `Opacity` sets how see-through that structure is; parts hidden in `Layers…` stay hidden and a part with its own opacity there keeps it.
- `Current` → `Contact`: pick a contact or net and the panel 'Current under a contact' shows the current density through it node by node - a profile across it (2-D sections) or a map over it (3-D devices) - with J max/mean and centre/edge. Emitter current crowding shows as a profile that rises towards the emitter edges.
- `I-V` BJT panel: β (left axis) and f_T (right axis, dotted; dash-dot for pins) against I_C for Gummel sweeps; the stars mark f_T,max.
- 3D Slicer `log10 G avalanche`: where impact ionisation makes electron-hole pairs (cm⁻³s⁻¹).
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
- Doping-dependent mobility (Masetti in Si, Caughey-Thomas/Arora form elsewhere) with velocity saturation.
- SRH with doping-dependent lifetimes, Auger and radiative recombination.
- Deep levels of defects: SRH through any level with its own capture cross sections, plus the charge it traps (see Defects).
- Oxide gates with fixed charge, interface sheet charges (GaN polarisation) and buried metal electrodes.
## Nanoscale MOS models (Solver page; MOS templates switch them on)
- Quantum confinement, MLDA (modified local-density approximation): at every semiconductor/insulator interface the density of states is cut by the factor 1 − exp(−(z/λ)²), λ = ħ/√(2 m kT), per valley: Si electrons in the Δ2 (m = 0.916, λ = 1.27 nm, 1/3 of the states) and Δ4 valleys (m = 0.19, λ = 2.79 nm, 2/3), holes heavy/light. The inversion layer moves ~1 nm off the oxide: V_th rises by 40-100 mV and the gate capacitance drops. The Bands tab draws Ec + Λn and Ev − Λp (the quantum potentials).
- Lombardi surface mobility (Si): 1/µ = 1/µ_bulk + D/µ_ac + D/µ_sr with the field normal to the current, D = exp(−d/10 nm). Reproduces the universal mobility curve (~400 cm²/Vs at 0.5 MV/cm, ~250 at 1 MV/cm for electrons).
- Gate tunnelling: Tsu-Esaki current with WKB transmission through the gate stack for electrons (conduction band) and holes (valence band), both directions (semiconductor → gate and gate → semiconductor). Direct tunnelling and Fowler-Nordheim are the same formula. Barriers from the electron affinities: SiO2 3.1 eV (tunnelling mass 0.42), HfO2 1.5 eV (0.18), Al2O3 2.7 eV (0.35), Si3N4 2.0 eV (0.5). The supply of an inversion or accumulation layer comes from its pressure on the interface, so the result does not depend on the mesh or on MLDA. The gate current appears in the Results table and the I-V semi-log panel; the 3D Slicer shows log10|J gate| at the interface nodes.
## Impact ionisation (Solver page; off by default, on in the avalanche templates)
- Electron-hole pairs are generated at G = (α_n|J_n| + α_p|J_p|)/q, with α = γ a exp(−(γ b/F)^m) (Chynoweth law), F the field along each carrier's current, and γ = tanh(ħω/2kT₀)/tanh(ħω/2kT) the optical-phonon temperature factor (ħω = 63 meV, T₀ = 300 K): breakdown voltages rise with temperature.
- Coefficients: Si van Overstraeten - de Man (holes: two field ranges, split at 4e5 V/cm); GaAs Bulman et al. (m = 1.6 / 1.75); 4H-SiC Hatakeyama et al. (along the c axis); GaN Ji et al. Germanium, InP and the alloys have none (no avalanche there).
- The generation is a source in both continuity equations (and in the balance of floating regions), updated every Gummel iteration.
- Near breakdown the multiplication grows without bound (M ~ 10-50) and the Gummel iteration, which solves the equations one after the other, stops converging; the I-V tab takes BV where 1/I, extrapolated from the last converged points, reaches zero. Past breakdown the current is held back only by the series resistance and the space charge of the avalanche carriers, and the solver may converge again on that high-current branch (a jump in the I-V curve).
- Checked in validate.py section 17: Si p+/n junctions break down where the ionisation integral of the simulated field reaches 1 (within 0.2 %, at 300 and 400 K), 3 % under Baliga's closed form.
## Quasi-static f_T (I-V tab, BJT panel)
- Every sweep point stores the electron and hole charge in the device. Along a Gummel sweep (V_CE held) f_T = g_m / (2π C_in), g_m = dI_C/dV_BE, C_in = dQ/dV_BE with Q the base's majority carrier (holes for an NPN): the charge-control estimate of where |h21| falls to 1. It contains the emitter and collector depletion charges and the stored minority charge (base, emitter, and the collector once the base pushes into it - the Kirk effect); it leaves out the non-quasi-static phase delay and every parasitic outside the simulated cell.
## Metal contacts (Solver page; on by default)
- Thermionic emission: at a Schottky contact or a contact with a metal, carriers are not pinned to their equilibrium densities; electrons and holes cross the interface with the Richardson velocity v = A**T²/(q·N_c,v) (A** = 112 A/cm²K² for n-Si, 32 for p-Si; per material) against the density in equilibrium with the metal. Drift-diffusion then carries them through the depletion region: together this is the thermionic-emission-diffusion theory of Crowell and Sze.
- Image-force lowering: the field at the interface lowers the barrier by Δφ = √(qF/4πε) (Schottky effect); the lowered barrier - flat up to the maximum x_m = Δφ/2F, −q/16πεx beyond - is part of the band the carriers move in, screened by the free carriers (Debye length) so that it acts in the depletion region only. Reverse currents therefore keep rising with voltage, as on real diodes.
- Tunnelling through the barrier (field and thermionic-field emission): Tsu-Esaki current with WKB transmission through the band along the contact normal (60 nm deep), injected where the carriers leave the barrier; it is what makes a metal on heavily doped silicon ohmic.
- Velocity saturation is tapered off for motion across the barrier region of a metal contact (carriers emitted over or tunnelling through it are hot); motion along the contact (a channel under a Schottky gate) keeps it.
- Off: Schottky contacts are the classic fixed-density barrier and Ohmic contacts with a metal stay ideal.
## Solver
- Newton-Poisson with multigrid-preconditioned CG.
- Continuity with BiCGSTAB and aggregation multigrid (on a prism mesh: strength-based pairwise aggregation of the mesh graph, aggregates of up to 8 nodes along the strongest couplings).
- A Gummel map on (ψ, φn, φp) accelerated by Anderson mixing.
- A floating-region balance step for regions tied to no contact.
- Adaptive bias continuation: voltage steps grow while they converge easily and shrink when they fail; inside a ramp each step starts from the extrapolation of the last two.
- Sweep predictor: when a solve continues the previous bias change in the same direction (the next point of a sweep, or the secant steps of a floating terminal), it starts from the secant extrapolation of the states before and after the previous solve instead of from the old state - 20-30 % fewer Gummel iterations on typical sweeps, about half on the IGBT; the results are unchanged (to ~1e-5).
- Terminal currents are exact in the discrete Kirchhoff sense and come with a round-off resolution floor. Currents below the floor are shown hollow or omitted.
## Validation
`validate.py` compares the core with closed-form theory (section 15, defects: see Defects; section 16: adaptive mesh refinement against uniformly refined meshes, and the sweep predictor against plain continuation): built-in voltage, depletion width, Shockley diode current, SiC blocking field, heterojunction offsets, MOS threshold, MOSFET linear current, the BJT Gummel integral and more. Section 14 checks the metal contacts: Schottky diodes against Crowell-Sze thermionic-emission-diffusion theory with image-force lowering (forward and reverse, within 0.3 %), the specific contact resistivity of a 0.65 eV barrier on 1e18-2e20 cm⁻³ n-Si (7 decades) against an independent Tsu-Esaki calculation, and aluminium being ohmic on n+ and rectifying on n− silicon. Section 13 checks the triangular-prism mesh: with every rectangular node kept it reproduces the rectangular solution exactly, the thinned mesh agrees within 1 % on MOSFETs, a GaN HEMT and the FinFET, and a round MOS capacitor (Si wire in an oxide shell) matches the radial Poisson-Boltzmann solution within a few mV where a staircase mesh misses by 20-80 mV. Sections 11-12 check the tunnelling kernel against direct integration, gate leakage against textbook values (1.2 nm SiO2: 1e2-1e3 A/cm² at 1 V; ~5 decades per nm), the Fowler-Nordheim slope against measured Si/SiO2 data, the MLDA threshold shift and inversion-layer centroid against quantum-mechanical results, and Lombardi against the universal mobility curve.
"""),
("Limitations", """
# Limitations and troubleshooting
## Not modelled
- impact ionisation is a local-field model (no dead space, no energy relaxation): very short high-field regions (below ~100 nm) break down at too low a voltage; the voltage-driven contacts make snapback and second breakdown out of reach; the runaway shows only when the leakage it multiplies is above the current resolution (in tiny cells with long lifetimes shorten τ0 on the Solver page - it raises the leakage, not BV); no coefficients for Ge, InP and the alloys;
- band-to-band tunnelling: Zener behaviour and GIDL are absent; it also adds to avalanche breakdown above ~5e16 cm⁻³, and the silicon ionisation coefficients (van Overstraeten & de Man) are fitted for 1.75–6e5 V/cm only - an abrupt Si junction at 1e17 cm⁻³ breaks down at 14 V here against about 11 V measured;
- incomplete ionisation of dopants: every donor and acceptor is ionised, so p-type 4H-SiC (Al) and p-GaN (Mg), where only a few per cent are ionised at 300 K, are far too conductive;
- Fermi-Dirac statistics: carriers follow Boltzmann statistics, so degenerate layers (n++/p++ contacts, GaAs above ~4e17 cm⁻³, HEMT channels) are approximate;
- negative differential mobility: GaAs and InP electrons saturate smoothly at v_sat (no transferred-electron velocity peak, no Gunn effect), so short GaAs/InP channels carry 2–3x too little current near the peak field;
- temperature dependence: silicon's saturation velocity follows Canali et al.; the other materials keep their 300 K v_sat, and the SiC/GaN Auger and ionisation coefficients have no temperature dependence;
- currents near avalanche: the multiplied current is converged to its resolution floor; below ~1e-16 A per cell it can still move by a few per cent between solves (the floor shown in the results says by how much);
- optical generation: solar cells show dark I-V only;
- quantisation beyond MLDA: no subband energies, and MLDA underestimates the dark space at very high fields; no quantum confinement in heterostructure 2DEGs;
- gate tunnelling: no image-force barrier lowering, no trap-assisted tunnelling, no valence-electron (EVB) tunnelling;
- thermionic emission at heterointerfaces;
- ballistic transport, strain and remote-phonon scattering in high-k stacks: nanoscale on-currents are drift-diffusion estimates;
- self-heating;
- transients and AC (f_T is the quasi-static charge-control estimate);
- defects: traps are in steady state (no transient trapping, current collapse or hysteresis), without field-enhanced emission (Poole-Frenkel, trap-assisted tunnelling); a metal particle carries no current through itself (a particle spanning a junction does not short it), and dislocations are not conducting pipes;
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
- File → Save device / Open device: the full editor state (regions, contacts, nets, sheets, defects, mesh, solver and sweep settings) as JSON.
- File → Export solution CSV: every node's x, y, z, φ, n, p, Ec, Ev, J and mobility.
- File → Export sweep CSV: the voltage, the floating-node voltage and the current and resolution floor of every terminal.
- File → Export plot: the visible tab as PNG, PDF or SVG.
- `validate.py`: physics regression battery (python3 validate.py; --templates also solves every template).
- `install.sh`: installs what is missing, builds the core and adds a launcher (Linux, WSL2, macOS).
- `semibuild.py`: builds the C core (python3 semibuild.py; --portable for any CPU of this kind). main.py runs it by itself when the core is missing or out of date.
- `semimesh.py`: the triangular-prism mesher; `semidefects.py`: the defect library and its mapping onto a mesh. Both must sit next to main.py.
"""),
("Keyboard shortcuts", "# Keyboard shortcuts\n"+"\n".join(
    f"  {k.replace('{MOD}',MOD):<14}{v}" for k,v in [
        ("F1","Help"),
        ("F5","Run (solve at the contact voltages)"),
        ("F6","Sweep"),
        ("Esc","Stop the running solve or sweep"),
        ("F11","Fill the screen / restore the window size"),
        ("{MOD}+E","Show / hide the editor panel"),
        ("{MOD}+O","Open device (JSON)"),
        ("{MOD}+S","Save device (JSON)"),
        ("{MOD} +/-/0","Text size larger / smaller / reset"),
        ("{MOD}+Q","Quit"),
        ("Mouse wheel","Scroll the panel under the pointer"),
        ("Shift+wheel","Scroll sideways"),
    ]+([("Cmd+R","Run (Mac keyboards need fn for F5)"),("Cmd+.","Stop")] if IS_MAC else []))+"\n"),
("Sources", _sources_text()),
]

def _sources_text():
    """docs/SOURCES.md - the literature and software SemiSim 3D is built on."""
    p=os.path.join(os.path.dirname(os.path.abspath(__file__)),"docs","SOURCES.md")
    try:
        with open(p,encoding="utf-8") as f: return f.read()
    except Exception:
        return ("# Sources\nThe list of sources is in docs/SOURCES.md in the SemiSim 3D folder "
                "(and on the project's GitHub page).")

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

def defect_section(d, H, Vv, N, pos, L3, tol=1e-6):
    """How a defect shows in the section plane N = pos (axes names H, Vv, N):
    (cut, kind, geom) - kind 'circle' geom (h, v, r) or 'rect' (h0, v0, h1, v1);
    cut: the plane cuts the defect (else geom is its projection on the plane).
    None: no outline (everywhere, or a slab parallel to the plane missing it)."""
    if d.shape=="everywhere": return None
    ih,iv,inn="xyz".index(H),"xyz".index(Vv),"xyz".index(N)
    lo,hi=d.extent(L3); c=np.array([d.x,d.y,d.z],float); r=0.5*d.a
    dn=pos-c[inn]
    if d.shape=="sphere":
        cut=abs(dn)<=r+tol
        return cut,"circle",(c[ih],c[iv],float(np.sqrt(max(r*r-dn*dn,0.))) if cut else r)
    if d.shape=="box":
        return lo[inn]-tol<=pos<=hi[inn]+tol,"rect",(lo[ih],lo[iv],hi[ih],hi[iv])
    ax="xyz".index(d.axis)
    if d.shape=="cylinder":
        if ax==inn:
            return lo[inn]-tol<=pos<=hi[inn]+tol,"circle",(c[ih],c[iv],r)
        cut=abs(dn)<=r+tol; w=float(np.sqrt(max(r*r-dn*dn,0.))) if cut else r
        if ax==ih: return cut,"rect",(lo[ih],c[iv]-w,hi[ih],c[iv]+w)
        return cut,"rect",(c[ih]-w,lo[iv],c[ih]+w,hi[iv])
    if ax==inn:                                          # slab parallel to the plane
        return (True,"rect",(lo[ih],lo[iv],hi[ih],hi[iv])) if abs(dn)<=r+tol else None
    cut=lo[inn]-tol<=pos<=hi[inn]+tol
    if ax==ih: return cut,"rect",(c[ih]-r,lo[iv],c[ih]+r,hi[iv])
    return cut,"rect",(lo[ih],c[iv]-r,hi[ih],c[iv]+r)

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
    if IS_MAC:
        # free + inactive + speculative + purgeable pages (what macOS hands out
        # without swapping); macOS also compresses memory, so never less than
        # half of the RAM
        try: half=os.sysconf("SC_PAGE_SIZE")*os.sysconf("SC_PHYS_PAGES")//2
        except Exception: half=0
        try:
            out=subprocess.run(["vm_stat"],capture_output=True,text=True,timeout=5).stdout
            import re
            ps=int(re.search(r"page size of (\d+)",out).group(1)); tot=0
            for key in ("Pages free","Pages inactive","Pages speculative","Pages purgeable"):
                m=re.search(key+r":\s+(\d+)",out)
                if m: tot+=int(m.group(1))
            if tot>0 or half>0: return max(tot*ps,half)
        except Exception: pass
        if half>0: return half
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
        try:                                    # window / taskbar / Dock icon
            ic=os.path.join(os.path.dirname(os.path.abspath(__file__)),"docs","icon.png")
            if os.path.isfile(ic): self._icon=tk.PhotoImage(file=ic); self.iconphoto(True,self._icon)
        except Exception: pass
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
        # adaptive mesh refinement (rectangular mesh): settings and the added lines (nm)
        self.amr_on=tk.BooleanVar(value=False); self.amr_passes=tk.IntVar(value=3)
        self.amr_budget=tk.DoubleVar(value=2.); self.amr_kappa=tk.DoubleVar(value=1.)
        self.amr_nodes={}
        # solver (meanings: CG its / BiCGSTAB its / Gummel its per bias step /
        # Gummel tolerance (V) / relative terminal-current tolerance)
        self.T=tk.DoubleVar(value=300.); self.max_p=tk.IntVar(value=2000)
        self.max_c=tk.IntVar(value=1000); self.max_g=tk.IntVar(value=500)
        self.tol_p=tk.DoubleVar(value=1e-8); self.tol_c=tk.DoubleVar(value=1e-6)
        self.m_bgn=tk.BooleanVar(value=True); self.m_tau=tk.BooleanVar(value=True)
        self.m_vsat=tk.BooleanVar(value=True); self.reuse_sol=tk.BooleanVar(value=True)
        # nanoscale MOS physics (set by each template, see its "models" key)
        self.m_qc=tk.BooleanVar(value=False); self.m_lb=tk.BooleanVar(value=False)
        self.m_tun=tk.BooleanVar(value=False)
        # metal contacts: thermionic emission (+ image force), barrier tunnelling
        self.m_te=tk.BooleanVar(value=True); self.m_sbt=tk.BooleanVar(value=True)
        self.m_ii=tk.BooleanVar(value=False)          # impact ionisation (avalanche)
        self.tau0=tk.StringVar(value="")              # background SRH lifetime (s); "" = material defaults
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
        self.show_labels=tk.BooleanVar(value=True)    # ... with region names where they fit
        self.cur_axis=tk.StringVar(value="X"); self.cur_p1=tk.IntVar(value=0); self.cur_p2=tk.IntVar(value=0)
        self.cur_show_arrows2d=tk.BooleanVar(value=True); self.cur_show_arrows3d=tk.BooleanVar(value=True)
        # arrow density: arrows across the in-plane map / along x in the 3-D field
        def _ci(v,d,lo,hi):
            try: return max(lo,min(hi,int(v)))
            except Exception: return d
        self.cur_n2d=tk.IntVar(value=_ci(self.cfg.get("arrows2d"),22,*ARROWS_2D))
        self.cur_n3d=tk.IntVar(value=_ci(self.cfg.get("arrows3d"),8,*ARROWS_3D))
        self.fi_axis=tk.StringVar(value="X"); self.fi_p1=tk.IntVar(value=0); self.fi_p2=tk.IntVar(value=0)
        self.fi_show_field=tk.BooleanVar(value=True); self.fi_show_mob=tk.BooleanVar(value=True)
        # which panels each plot tab shows (Panels ▾ menus); the Structure 3-D
        # view and the Field tab's E / µ switches are the same variables
        self.pv={"td":{"3d":self.sec_3d},"tf":{"E":self.fi_show_field,"mu":self.fi_show_mob}}
        # Structure tab layers (Layers… window): hidden parts (by object id),
        # opacities, insulators in 3-D, defects, apply to the cross-section too
        self.lay_hide=set(); self.lay_def=tk.BooleanVar(value=True); self.lay_ins3d=tk.BooleanVar(value=False)
        self.lay_ra=tk.DoubleVar(value=0.35); self.lay_ca=tk.DoubleVar(value=0.85); self.lay_sec=tk.BooleanVar(value=False)
        self.lay_alpha={}                                 # id(part) -> its own opacity (Layers… window)
        self.lay_sa=tk.DoubleVar(value=0.35)              # the opacity slider of the selected parts
        self._lay_top=None
        # Current tab: the structure under the current arrows
        self.cur_bg=tk.StringVar(value="|J|")             # in-plane map: |J| / structure / |J| + structure
        self.cur_st3=tk.BooleanVar(value=False)           # structure in the 3-D current field
        self.cur_sa=tk.DoubleVar(value=0.25)              # its opacity
        # current under a contact (Current tab), pinned sweeps, A/B designs
        self.cur_prof=tk.StringVar(value="(none)")
        self.pins=[]; self._pin_top=None; self._ab_top=None; self.ab={"A":None,"B":None}
        self._ptab_top=None; self.ptab_ref=tk.StringVar(value="")       # the parameter table (I-V tab ▸ Table…)
        self.ptab_mode=tk.StringVar(value="Values and differences"); self.ptab_all=tk.BooleanVar(value=False)
        self.ab_keep=tk.BooleanVar(value=True)
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
        self.cb_metal=tk.StringVar(value=METAL_NONE)
        # defects (semidefects.py): the list, the defect form, random scatter
        self.defects=[]
        self.df_cat=tk.StringVar(value=SD.CATEGORIES[0]); self.df_kind=tk.StringVar(value="Void")
        self.df_shape=tk.StringVar(value="sphere"); self.df_axis=tk.StringVar(value="y")
        self.df_x=tk.StringVar(value="0"); self.df_y=tk.StringVar(value="0"); self.df_z=tk.StringVar(value="0")
        self.df_a=tk.StringVar(value="100"); self.df_b=tk.StringVar(value="0"); self.df_c=tk.StringVar(value="0")
        self.df_dens=tk.StringVar(value="0"); self.df_dit=tk.StringVar(value="0"); self.df_q=tk.StringVar(value="0")
        self.df_metal=tk.StringVar(value="NiSi"); self.df_host=tk.StringVar(value="(any)")
        self.df_lt=tk.StringVar(value="A"); self.df_lr=tk.StringVar(value="Ei"); self.df_lE=tk.StringVar(value="0")
        self.df_lsn=tk.StringVar(value="1e-15"); self.df_lsp=tk.StringVar(value="1e-15")
        self.df_lbl=tk.StringVar(value="")
        self.sc_mode=tk.StringVar(value="count"); self.sc_n=tk.StringVar(value="20")
        self.sc_seed=tk.StringVar(value="1"); self.sc_spread=tk.StringVar(value="0")
        self.sc_box=[tk.StringVar(value="0") for _ in range(6)]        # x0 x1 y0 y1 z0 z1 (nm)
        self.sc_grp=tk.StringVar(value="")
        # text size
        sc=self.cfg.get("text_scale",1.0)
        self.txt_scale=tk.DoubleVar(value=sc if sc in self.TEXT_SCALES else 1.0)
        # defaults for widgets that do not set their own colours (the
        # matplotlib toolbars pick light icons on a dark button background)
        opts=[("*Label.background","#020609"),("*Label.foreground","#9fb3c8")]
        if not _aqua(self):
            # (macOS draws these natively in light colours: light icons would vanish)
            opts+=[("*Button.background","#0b1422"),("*Button.foreground","#c9d8ec"),
                   ("*Button.activeBackground","#182436"),("*Checkbutton.background","#0b1422"),
                   ("*Checkbutton.foreground","#c9d8ec"),("*Checkbutton.selectColor","#182436")]
        for k,v in opts: self.option_add(k,v,"widgetDefault")
        self._win_init(); _mac_no_app_nap()
        self._fonts(); self._style(); self._build_ui(); self._bind_keys()
        self._apply_text_scale(redraw=False)
        self.protocol("WM_DELETE_WINDOW",self._quit)
        last=self.cfg.get("template")
        self._load(last if last in TEMPLATES else "PN Diode")
        if not self.cfg.get("drawer",True): self.after(50,self._toggle_drawer)
        if self.cfg.get("filled"): self.after(400,self._toggle_fill)    # it was closed full-size
        if IS_MAC:             # started from a terminal: come to the front
            self.lift(); self.attributes("-topmost",True)
            self.after(600,lambda:self.attributes("-topmost",False)); self.focus_force()

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
            # a maximised, filled or full-screen window is saved as the size it returns to
            big=self._fill_on or self._zoomed() or self._fullscreen()
            geo=(self._normal_geo or self._restore_geo or self.geometry()) if big else self.geometry()
            d=dict(geometry=geo,text_scale=self.txt_scale.get(),
                   template=self.tname,drawer=self._drawer_on,drawer_w=w or self._drawer_w,
                   workarea=self._workarea,wsl_max_fix=bool(self.wsl_fix.get()),filled=bool(big),
                   arrows2d=self._arrow_n()[0],arrows3d=self._arrow_n()[1])
            with open(SETTINGS_FILE,"w") as f: json.dump(d,f)
        except Exception: pass
    def _quit(self):
        if self.busy and self.job is not None: self.job.request_stop()
        self._save_settings(); self.destroy()

    # ─── window size: fill the screen (F11) ────────────────────
    # WSLg (the window system of WSL2) mis-routes the mouse in MAXIMISED X11
    # windows: after the maximise button, clicks miss the buttons while a
    # drag on a plot still rotates it.  A normal window as large as the
    # screen does not have that problem, so on WSL F11 makes one, and a
    # maximise (button, title double-click, Win+Up) is turned into one: the
    # window is un-maximised and given the rectangle the maximise chose.
    def _win_init(self):
        self._fill_on=False; self._restore_geo=None; self._unzoom_job=None
        self._normal_geo=None; self._last_cmd={}
        wa=self.cfg.get("workarea")
        self._workarea=wa if isinstance(wa,str) and wa.count("+")==2 else None
        self.wsl_fix=tk.BooleanVar(value=IS_WSL and bool(self.cfg.get("wsl_max_fix",True)))
        self.bind("<Configure>",self._win_configure,add="+")

    def _zoomed(self):
        """Is the window maximised by the window manager?"""
        try:
            if _x11(self):
                return bool(int(self.wm_attributes("-zoomed")))
            return str(self.wm_state())=="zoomed"
        except Exception: return False

    def _fullscreen(self):
        try: return bool(int(self.wm_attributes("-fullscreen")))
        except Exception: return False

    def _win_configure(self,e):
        if e.widget is not self: return
        if self._zoomed():
            if IS_WSL and self.wsl_fix.get() and self._unzoom_job is None:
                self._unzoom_job=self.after(250,self._wsl_unzoom)
        elif not self._fill_on:
            # remember the normal size to return to, once it has held for half
            # a second (a maximise passes through intermediate sizes first)
            self._debounce("_normal_job",self._commit_normal,500)

    def _commit_normal(self):
        if self._fill_on or self._zoomed() or self._fullscreen(): return
        if (self.winfo_width()<self.winfo_screenwidth()-8 or
                self.winfo_height()<self.winfo_screenheight()-160):    # not a screen-sized one
            self._normal_geo=self.geometry()

    def _wsl_unzoom(self):
        self._unzoom_job=None
        if not self._zoomed(): return
        now=time.monotonic(); hist=[t for t in getattr(self,"_unzoom_t",[]) if now-t<6.]+[now]
        self._unzoom_t=hist
        if len(hist)>3:          # the window system insists on maximising: leave it alone
            self.wsl_fix.set(False)
            self.sv.set("Maximised by the window system - if clicks stop working, press F11 twice")
            return
        self.update_idletasks()
        geo=self.geometry()                     # the rectangle the maximise chose = the work area
        if not self._fill_on: self._restore_geo=self._normal_geo
        self._workarea=geo
        try: self.wm_attributes("-zoomed",0)
        except Exception: return
        self._fill_on=True
        self.after(150,lambda:self._fill_apply(geo))
        self.sv.set("Window filled to the screen (F11 restores) - WSLg mis-routes clicks in maximised windows")

    def _fill_apply(self,geo,tries=3):
        """Set geometry and check it stuck (the window manager may still be
        restoring the old size when the first request arrives)."""
        try: self.geometry(geo)
        except Exception: return
        def check():
            import re
            m=re.match(r"(\d+)x(\d+)",geo)
            if m and tries>1 and self._fill_on and not self._zoomed():
                w,h=int(m.group(1)),int(m.group(2))
                if abs(self.winfo_width()-w)>4 or abs(self.winfo_height()-h)>4:
                    self._fill_apply(geo,tries-1)
        self.after(350,check)

    def _screen_geo(self):
        """The screen area for F11: learnt from a maximise when there was one,
        else the screen minus the Windows taskbar and the title bar.  Several
        monitors side by side form one wide screen: use the one under the window."""
        if self._workarea: return self._workarea
        sw,sh=self.winfo_screenwidth(),self.winfo_screenheight()
        k=max(1,int(round(sw/(sh*16/9.)))) if sw>2*sh else 1
        mw=sw//k; cx=self.winfo_x()+self.winfo_width()//2
        i=min(k-1,max(0,cx//max(mw,1)))
        top=0; bottom=48+34 if IS_WSL else 60
        return f"{mw}x{max(400,sh-bottom)}+{i*mw}+{top}"

    def _toggle_fill(self):
        if self._fill_on or self._zoomed():
            if self._zoomed():
                try: (self.wm_attributes("-zoomed",0) if _x11(self)
                      else self.wm_state("normal"))
                except Exception: pass
            self._fill_on=False
            g=self._normal_geo or self._restore_geo      # the last normal size
            if g: self.after(150,lambda:self.geometry(g))
            return
        self._restore_geo=self._normal_geo or self.geometry()
        if IS_WSL:
            self._fill_on=True; self._fill_apply(self._screen_geo())
            return
        try:
            if _x11(self): self.wm_attributes("-zoomed",1)
            else: self.wm_state("zoomed")
        except Exception:
            self._fill_on=True; self._fill_apply(self._screen_geo())

    def _reset_window(self):
        """Back to the default size in the middle of the screen."""
        if self._zoomed():
            try: (self.wm_attributes("-zoomed",0) if _x11(self)
                  else self.wm_state("normal"))
            except Exception: pass
        self._fill_on=False; self._workarea=None
        sw,sh=self.winfo_screenwidth(),self.winfo_screenheight()
        W=min(1680,max(800,sw-40)); H=min(1020,max(560,sh-90))
        g=f"{W}x{H}+{max(0,(sw-W)//2)}+{max(0,(sh-H)//4)}"
        self.after(150,lambda:self.geometry(g))

    def _once(self,key,fn):
        """fn behind a 0.35 s guard: on macOS a menu shortcut can arrive both
        as the menu's key equivalent and as a key binding."""
        def w(*_):
            t=time.monotonic()
            if t-self._last_cmd.get(key,-1.)<0.35: return "break"
            self._last_cmd[key]=t; fn(); return "break"
        return w

    # ─── fonts & style ─────────────────────────────────────────
    def _fonts(self):
        ui=_first_family(self,["Segoe UI","Ubuntu","DejaVu Sans","Noto Sans","Liberation Sans",
                               "Helvetica Neue","Arial","Helvetica"],"TkDefaultFont")
        mono=_first_family(self,["Consolas","Cascadia Mono","DejaVu Sans Mono","Ubuntu Mono",
                                 "Liberation Mono","Noto Sans Mono","Menlo","Courier New"],"TkFixedFont")
        # Tk on macOS counts 72 points per inch of the screen (96 elsewhere):
        # the same point size would come out 25 % smaller there
        k=1.25 if _aqua(self) else 1.0
        mk=lambda fam,sz,**kw: tkfont.Font(self,family=fam,size=int(round(sz*k)),**kw)
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
        cls=FlatButton if FLAT_BUTTONS else tk.Button
        return cls(p,text=t,command=cmd,bg=bg or C["accent"],fg=fg,activebackground=C["border"],
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
        cv.mpl_connect("draw_event",self._clear_3d_on_draw)
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
        aqua=_aqua(self)
        if aqua:            # macOS application menu (SemiSim 3D ▸ About, Quit)
            am=tk.Menu(mb,name="apple",tearoff=0)
            am.add_command(label="About SemiSim 3D",command=self._about); am.add_separator()
            mb.add_cascade(menu=am)
            try: self.createcommand("tk::mac::Quit",self._quit)
            except Exception: pass
        A=self._acc_cmd
        fm=tk.Menu(mb,**mkw)
        fm.add_command(label="Open device…",accelerator=f"{MOD}+O",command=A("open",self._open_device))
        fm.add_command(label="Save device as…",accelerator=f"{MOD}+S",command=A("save",self._save_device))
        fm.add_separator()
        fm.add_command(label="Export solution CSV…",command=self._xcsv)
        fm.add_command(label="Export sweep CSV…",command=self._xiv_csv)
        fm.add_command(label="Export visible plot…",command=self._xplt)
        if not aqua:        # (macOS: Quit is in the application menu)
            fm.add_separator(); fm.add_command(label="Quit",accelerator=f"{MOD}+Q",command=self._quit)
        mb.add_cascade(label="File",menu=fm)
        dm=tk.Menu(mb,**mkw); self._fill_template_menu(dm,mkw)
        mb.add_cascade(label="Devices",menu=dm)
        rm=tk.Menu(mb,**mkw)
        rm.add_command(label="▶ Run",accelerator="F5",command=A("run",self._run))
        rm.add_command(label="↔ Sweep",accelerator="F6",command=A("sweep",self._sweep))
        rm.add_command(label="■ Stop",accelerator="Esc",command=self._stop)
        rm.add_separator()
        rm.add_checkbutton(label="Continue from previous solution",variable=self.reuse_sol)
        rm.add_command(label="Forget solution (next run starts from equilibrium)",command=self._forget)
        mb.add_cascade(label="Run",menu=rm)
        vm=tk.Menu(mb,**mkw)
        vm.add_command(label="Show / hide editor panel",accelerator=f"{MOD}+E",command=A("drawer",self._toggle_drawer))
        vm.add_command(label="Fill the screen / restore",accelerator="F11",command=A("fill",self._toggle_fill))
        vm.add_command(label="Reset window size",command=self._reset_window)
        if IS_WSL:
            vm.add_checkbutton(label="WSL: turn maximise into fill screen (fixes dead clicks)",
                               variable=self.wsl_fix)
        vm.add_separator()
        ts=tk.Menu(vm,**mkw)
        for sc in self.TEXT_SCALES:
            ts.add_radiobutton(label=f"{int(sc*100)} %",variable=self.txt_scale,value=sc,
                               command=self._apply_text_scale)
        vm.add_cascade(label="Text size",menu=ts)
        vm.add_command(label="Larger text",accelerator=f"{MOD}++",command=A("txt+",lambda:self._text_step(1)))
        vm.add_command(label="Smaller text",accelerator=f"{MOD}+-",command=A("txt-",lambda:self._text_step(-1)))
        vm.add_separator()
        for i,(nm,_) in enumerate(self._TABS):
            vm.add_command(label=f"{nm.strip()} tab",command=lambda i=i:self.nb.select(i))
        mb.add_cascade(label="View",menu=vm)
        hm=tk.Menu(mb,**mkw)
        hm.add_command(label="Help contents",accelerator="F1",command=self._help)
        hm.add_command(label="Keyboard shortcuts",command=lambda:self._help("Keyboard shortcuts"))
        hm.add_command(label="Template library…",command=self._browse_templates)
        hm.add_command(label="Status (test version) and terms",command=lambda:self._help("Status and terms"))
        hm.add_command(label="Sources",command=lambda:self._help("Sources"))
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
        Tip(b,f"Show / hide the editor panel ({MOD}+E)",self.F["small"])
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

    def _acc_cmd(self,key,fn):
        """Menu command of an entry that has a keyboard shortcut (see _once)."""
        if not _aqua(self): return fn
        w=self._once(key,fn); return lambda:w()

    def _bind_keys(self):
        aqua=_aqua(self)
        K=(lambda key,fn:self._once(key,fn)) if aqua else (lambda key,fn:(lambda e:(fn(),"break")[1]))
        self.bind_all("<F1>",lambda e:self._help())
        self.bind_all("<F5>",K("run",self._run))
        self.bind_all("<F6>",K("sweep",self._sweep))
        self.bind_all("<F11>",K("fill",self._toggle_fill))
        self.bind_all("<Escape>",lambda e:self._stop())
        def B(seq,fn):
            try: self.bind_all(seq,fn)
            except tk.TclError: pass        # (a Command binding outside macOS: tests)
        keys=[("e","drawer",self._toggle_drawer),("o","open",self._open_device),("s","save",self._save_device)]
        mods=["Control"]+(["Command"] if aqua else [])
        for mod in mods:
            for k,nm,fn in keys:
                for kk in (k,k.upper()): B(f"<{mod}-{kk}>",K(nm,fn))
            for seq in ("plus","equal","KP_Add"): B(f"<{mod}-{seq}>",K("txt+",lambda:self._text_step(1)))
            for seq in ("minus","KP_Subtract"): B(f"<{mod}-{seq}>",K("txt-",lambda:self._text_step(-1)))
            B(f"<{mod}-0>",K("txt0",lambda:self._text_step(None)))
        if aqua:            # Mac keyboards need fn for F5: Cmd+R runs, Cmd+. stops (Cmd+Q: application menu)
            B("<Command-r>",K("run",self._run))
            B("<Command-period>",lambda e:self._stop())
        else:
            B("<Control-q>",K("quit",self._quit))
        for seq in ("<MouseWheel>","<Shift-MouseWheel>","<Button-4>","<Button-5>",
                    "<Shift-Button-4>","<Shift-Button-5>"):
            self.bind_all(seq,self._wheel,add="+")
        try:                # Tk 9: trackpads send TouchpadScroll, not MouseWheel
            self.bind_all("<TouchpadScroll>",self._touch_scroll,add="+")
        except tk.TclError: pass

    def _wheel_target(self,e):
        """The ScrollFrame canvas under the pointer, or None (widgets that
        scroll by themselves - tables, text, plots, sliders - are left alone)."""
        try: x=self.winfo_containing(e.x_root,e.y_root)
        except Exception: return None
        while x is not None:
            if getattr(x,"wheel_native",False): return None
            if isinstance(x,(tk.Text,tk.Listbox,ttk.Treeview,ttk.Scale,tk.Scale)): return None
            if getattr(x,"wheel_target",False): return x
            x=getattr(x,"master",None)
        return None

    def _wheel(self,e):
        """Scroll the ScrollFrame under the pointer (Windows/macOS MouseWheel
        and X11 buttons 4/5; Shift = sideways)."""
        x=self._wheel_target(e)
        if x is None: return
        if getattr(e,"num",None)==4: step=-1
        elif getattr(e,"num",None)==5: step=1
        else:
            d=getattr(e,"delta",0)
            if not d: return
            step=-max(1,abs(d)//120) if d>0 else max(1,abs(d)//120)
        horiz=bool(e.state & 0x0001)
        (x.xview_scroll if horiz else x.yview_scroll)(step,"units"); return "break"

    def _touch_scroll(self,e):
        """Tk 9 trackpad scrolling (pixel deltas; macOS)."""
        x=self._wheel_target(e)
        if x is None: return
        try: dx,dy=(int(float(v)) for v in self.tk.splitlist(self.tk.call("tk::PreciseScrollDeltas",e.delta)))
        except Exception: return
        acc=self.__dict__.setdefault("_tacc",[0.,0.])
        acc[0]-=dx; acc[1]-=dy
        for i,view in ((0,x.xview_scroll),(1,x.yview_scroll)):
            n=int(acc[i]/24.)                   # the pages scroll in 24-px units
            if n: view(n,"units"); acc[i]-=24.*n
        return "break"

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
        for key,title in [("dev","Device"),("reg","Regions"),("con","Contacts"),("def","Defects"),
                          ("sol","Solver"),("swp","Sweep")]:
            sf=ScrollFrame(nb,C["panel"]); nb.add(sf,text=title); pages[key]=sf.inner
        self._page_device(pages["dev"]); self._page_regions(pages["reg"])
        self._page_contacts(pages["con"]); self._page_defects(pages["def"])
        self._page_solver(pages["sol"]); self._page_sweep(pages["swp"])

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
        fr=self._sec(p,"Adaptive mesh refinement")
        fm=tk.Frame(fr,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        ttk.Checkbutton(fm,text="Refine the mesh to the solution on ▶ Run (rectangular mesh)",
                        variable=self.amr_on).grid(row=0,column=0,columnspan=6,sticky="w",padx=4,pady=2)
        self._L(fm,"Passes",1,0); self._E(fm,self.amr_passes,1,1,5)
        self._L(fm,"Max nodes ×",1,2); self._E(fm,self.amr_budget,1,3,5)
        self._L(fm,"Threshold (kT/q)",2,0); self._E(fm,self.amr_kappa,2,1,5)
        self._btn(fm,"Clear refinement",self._amr_clear,bg=C["border"]).grid(row=2,column=2,columnspan=2,padx=2,sticky="w")
        self.amr_lbl=tk.Label(fr,text="no refinement",bg=C["panel"],fg=C["green"],font=self.F["small"],
                              anchor="w",justify="left")
        self.amr_lbl.pack(fill=tk.X,padx=8,pady=(0,2))
        self._hint(fr,"After the solve, every node is checked along each axis: where the potential, or the "
                      "quasi-Fermi potential of a carrier that matters there, bends away from the straight line "
                      "through its neighbours by more than the threshold (in kT/q = 26 mV), the mesh does not "
                      "resolve the solution - an unresolved space-charge region, an inversion layer, current "
                      "crowding. Those intervals are split (a new mesh line through the device) and the device is "
                      "solved again, starting from the interpolated solution; up to 'Passes' times, until "
                      "nothing is under-resolved or the mesh reaches 'Max nodes ×' the unrefined node count. "
                      "Uniform fields and resistors need no refinement and get none. The added lines stay for "
                      "sweeps and later runs (refine at the bias that matters most - e.g. the highest current - "
                      "then sweep) until 'Clear refinement' or a new template; they are saved with the device.")
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
        cols=("#","Label","Net","Face","BC","Metal","V","φm","EOT (stack)","Nox","Range (nm)")
        tf,self.ctree=self._tree(fr,cols,[30,150,60,56,60,70,56,46,90,52,260],height=8)
        tf.pack(fill=tk.BOTH,expand=True,padx=4,pady=4)
        self.ctree.bind("<<TreeviewSelect>>",lambda _:self._loadcon())
        fr2=self._sec(p,"Edit contact")
        fm=tk.Frame(fr2,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"Face",0,0); self._CB(fm,self.cb_face,Con.FNAME,0,1,w=8)
        self._L(fm,"BC",0,2); self._CB(fm,self.cb_bc,Con.BCNAME,0,3,w=9)
        self._L(fm,"V (V)",1,0); self._E(fm,self.cb_V,1,1,9)
        self._L(fm,"φm (eV)",1,2); self.cb_pm_e=self._E(fm,self.cb_pm,1,3,7)
        self._L(fm,"Metal",2,0); self._CB(fm,self.cb_metal,[METAL_NONE]+list(METALS),2,1,w=17,columnspan=3)
        self.metal_lbl=tk.Label(fm,text="",bg=C["panel"],fg=C["green"],font=self.F["small"],justify="left",anchor="w")
        self.metal_lbl.grid(row=3,column=0,columnspan=4,sticky="w",padx=4)
        self._L(fm,"Net",4,0); self._E(fm,self.cb_net,4,1,9)
        self._L(fm,"(same net = wired together)",4,2,columnspan=2)
        self.con_axis_lbl=tk.Label(fm,text="",bg=C["panel"],fg=C["yellow"],font=self.F["uib"])
        self.con_axis_lbl.grid(row=5,column=0,columnspan=4,sticky="w",padx=4,pady=(4,0))
        self.cb_face.trace_add("write",self._upd_con_axes)
        for v in (self.cb_metal,self.cb_bc,self.cb_face,self.cb_i0,self.cb_i1,self.cb_j0,self.cb_j1):
            v.trace_add("write",lambda *_:self._debounce("_metal_after",self._upd_metal,120))
        self.ci_lbl=self._L(fm,"i0 / i1 %",6,0); self._E(fm,self.cb_i0,6,1,8); self._E(fm,self.cb_i1,6,2,8)
        self.cj_lbl=self._L(fm,"j0 / j1 %",7,0); self._E(fm,self.cb_j0,7,1,8); self._E(fm,self.cb_j1,7,2,8)
        self.ck_lbl=self._L(fm,"k0 / k1 %",8,0); self._E(fm,self.cb_k0,8,1,8); self._E(fm,self.cb_k1,8,2,8)
        self._L(fm,"Gate tox (nm)",9,0); self._E(fm,self.cb_tox,9,1,8)
        self._L(fm,"Nox (q/cm²)",9,2); self._E(fm,self.cb_nox,9,3,9)
        self._L(fm,"Gate stack",10,0); self._E(fm,self.cb_stack,10,1,26,columnspan=3)
        self._L(fm,"Label",11,0); self._E(fm,self.cb_lbl,11,1,26,columnspan=3)
        bf=tk.Frame(fr2,bg=C["panel"]); bf.pack(fill=tk.X,padx=4,pady=4)
        items=[("+ Add",self._addcon,"#0a8f5a"),("✎ Update",self._editcon,C["accent"]),
               ("⧉ Duplicate",self._dupcon,"#0e7490"),("✕ Delete",self._delcon,"#a3243b"),
               ("Full face",self._full_face,C["border"]),("Clear all",self._clrcon,"#444")]
        for k,(t2,cmd,bg) in enumerate(items):
            self._btn(bf,t2,cmd,bg=bg).grid(row=k//3,column=k%3,padx=2,pady=2,sticky="ew")
        for c in range(3): bf.columnconfigure(c,weight=1)
        self._hint(fr2,"Ranges are % of the face's two in-plane axes (Box: % of X, Y and Z). "
                       "Ohmic = charge-neutral contact; Schottky uses φm; Metal: the contact is made of that "
                       "metal - its measured Schottky barrier on the semiconductor underneath (or the Cowley-Sze "
                       "pinning model) replaces φm, and thermionic emission plus tunnelling through the barrier "
                       "decide whether an 'Ohmic' contact really is ohmic (heavy doping) or rectifies (Solver "
                       "page, metal-contact physics); a gate takes the metal's work function. Gate on a semiconductor face = "
                       "oxide of thickness tox (EOT) with fixed charge Nox; Gate on oxide nodes or as a Box = "
                       "metal electrode. Gate stack (face gates, optional): layers from the semiconductor "
                       "outwards in nm, e.g. 'SiO2 0.5, HfO2 2.5' - it replaces tox (EOT computed from the "
                       "permittivities) and is what gate tunnelling sees; empty = SiO2 of thickness tox. "
                       "Contacts sharing a Net always have the same voltage - editing one member's V "
                       "updates the whole net.")
        self._upd_con_axes()

    def _page_defects(self,p):
        fr=self._sec(p,"Defects (applied in order)")
        cols=("#","Label","Kind","Shape (nm)","Centre x, y, z (nm)","Strength","Host")
        tf,self.dtree=self._tree(fr,cols,[34,190,190,110,150,130,60],height=8)
        tf.pack(fill=tk.BOTH,expand=True,padx=4,pady=4)
        self.dtree.bind("<<TreeviewSelect>>",lambda _:self._loaddef())
        self.def_sum=tk.Label(fr,text="no defects",bg=C["panel"],fg=C["green"],font=self.F["small"],
                              justify="left",anchor="w",wraplength=self._hint_wrap())
        self.def_sum.pack(fill=tk.X,padx=8,pady=(0,4)); self._hints=getattr(self,"_hints",[])+[self.def_sum]
        fr2=self._sec(p,"Edit defect")
        fm=tk.Frame(fr2,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"Category",0,0)
        cb=self._CB(fm,self.df_cat,SD.CATEGORIES,0,1,w=30,columnspan=4)
        cb.bind("<<ComboboxSelected>>",lambda e:self._def_cat_changed())
        self._L(fm,"Defect",1,0)
        self.df_kind_cb=self._CB(fm,self.df_kind,[],1,1,w=30,columnspan=4)
        self.df_kind_cb.bind("<<ComboboxSelected>>",lambda e:self._def_kind_changed())
        self.df_desc=tk.Label(fm,text="",bg=C["panel"],fg=C["text"],font=self.F["small"],justify="left",
                              anchor="w",wraplength=self._hint_wrap())
        self.df_desc.grid(row=2,column=0,columnspan=5,sticky="w",padx=4,pady=(0,4))
        self._hints.append(self.df_desc)
        self._L(fm,"Shape",3,0)
        self.df_shape_cb=self._CB(fm,self.df_shape,SD.SHAPES,3,1,w=10)
        self.df_shape_cb.bind("<<ComboboxSelected>>",lambda e:self._def_fields())
        self._L(fm,"Axis",3,2); self.df_axis_cb=self._CB(fm,self.df_axis,["x","y","z"],3,3,w=3)
        self._L(fm,"Centre x / y / z",4,0)
        self.df_ce=[self._E(fm,v,4,1+k,8) for k,v in enumerate((self.df_x,self.df_y,self.df_z))]
        self._btn(fm,"centre",self._def_centre,bg=C["border"]).grid(row=4,column=4,padx=2,sticky="w")
        self.df_sz_lbl=self._L(fm,"Size (nm)",5,0)
        self.df_se=[self._E(fm,v,5,1+k,8) for k,v in enumerate((self.df_a,self.df_b,self.df_c))]
        self.df_sz_txt=tk.Label(fm,text="",bg=C["panel"],fg=C["sub"],font=self.F["small"],anchor="w")
        self.df_sz_txt.grid(row=6,column=1,columnspan=4,sticky="w",padx=2)
        self._L(fm,"Density",7,0); self.df_dens_e=self._E(fm,self.df_dens,7,1,10)
        self.df_dens_u=self._L(fm,"",7,2,columnspan=3)
        self._L(fm,"D_it (cm⁻²eV⁻¹)",8,0); self.df_dit_e=self._E(fm,self.df_dit,8,1,10)
        self._L(fm,"Charge",8,2); self.df_q_e=self._E(fm,self.df_q,8,3,9)
        self.df_q_u=self._L(fm,"",8,4)
        self._L(fm,"Metal",9,0); self.df_metal_cb=self._CB(fm,self.df_metal,[m for m in METALS if "poly" not in m],9,1,w=10)
        self._L(fm,"Host",9,2)
        self._CB(fm,self.df_host,["(any)"]+sorted(SEMIS,key=lambda m:MATS[m]),9,3,w=8,columnspan=2)
        self._L(fm,"Level",10,0)
        lf=tk.Frame(fm,bg=C["panel"]); lf.grid(row=10,column=1,columnspan=4,sticky="w")
        self.df_lv_w=[ttk.Combobox(lf,textvariable=self.df_lt,values=["A","D","N"],state="readonly",width=3,font=self.F["ui"]),
                      ttk.Combobox(lf,textvariable=self.df_lr,values=["Ec","Ev","Ei"],state="readonly",width=4,font=self.F["ui"])]
        for w in self.df_lv_w: w.pack(side=tk.LEFT,padx=2)
        e=tk.Entry(lf,textvariable=self.df_lE,bg="#07111e",fg=C["text"],insertbackground=C["text"],relief="flat",
                   font=self.F["mono"],width=7,highlightbackground=C["border"],highlightcolor=C["accent"],highlightthickness=1)
        e.pack(side=tk.LEFT,padx=2); self.df_lv_w.append(e)
        tk.Label(lf,text="eV",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT)
        self._L(fm,"σn / σp (cm²)",11,0)
        self.df_lv_w+= [self._E(fm,self.df_lsn,11,1,9),self._E(fm,self.df_lsp,11,2,9)]
        self.df_lvl_txt=tk.Label(fm,text="",bg=C["panel"],fg=C["green"],font=self.F["small"],justify="left",
                                 anchor="w",wraplength=self._hint_wrap())
        self.df_lvl_txt.grid(row=12,column=0,columnspan=5,sticky="w",padx=4); self._hints.append(self.df_lvl_txt)
        self._L(fm,"Label",13,0); self._E(fm,self.df_lbl,13,1,30,columnspan=4)
        bf=tk.Frame(fr2,bg=C["panel"]); bf.pack(fill=tk.X,padx=4,pady=4)
        items=[("+ Add",self._adddef,"#0a8f5a"),("✎ Update",self._editdef,C["accent"]),
               ("⧉ Duplicate",self._dupdef,"#0e7490"),("✕ Delete",self._deldef,"#a3243b"),
               ("▲ Up",self._defup,C["border"]),("▼ Down",self._defdn,C["border"]),
               ("Clear all",self._clrdef,"#444")]
        for k,(t2,cmd,bg) in enumerate(items):
            self._btn(bf,t2,cmd,bg=bg).grid(row=k//4,column=k%4,padx=2,pady=2,sticky="ew")
        for c in range(4): bf.columnconfigure(c,weight=1)
        self._hint(fr2,"Pick a category and a defect: the library values (literature, editable) fill the form. "
                       "Geometry in nm (y = depth); the centre is ignored by 'everywhere'. Size: sphere = "
                       "diameter; box = dx, dy, dz; cylinder = diameter and length along its axis (0 = through the "
                       "device); slab = thickness along its normal axis and the two widths (0 = across the device). "
                       "Deep levels: acceptor (−) when filled, donor (+) when empty, neutral = recombination only; "
                       "energy below Ec, above Ev or from midgap (Ei); editing the level replaces the first level "
                       "of the defect. Host = the semiconductor the defect lives in. Inclusions (voids, oxide "
                       "and metal particles) replace the material in their shape; a metal particle pins the Fermi "
                       "level at its Schottky barrier. Later rows override earlier ones where inclusions overlap.")
        fr3=self._sec(p,"Random scatter")
        fm=tk.Frame(fr3,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        ttk.Radiobutton(fm,text="Number",variable=self.sc_mode,value="count",
                        command=self._def_sc_unit).grid(row=0,column=0,sticky="w",padx=4)
        ttk.Radiobutton(fm,text="Density",variable=self.sc_mode,value="density",
                        command=self._def_sc_unit).grid(row=0,column=1,sticky="w")
        self._E(fm,self.sc_n,0,2,9); self.sc_unit=self._L(fm,"",0,3,columnspan=2)
        self._L(fm,"Seed",1,0); self._E(fm,self.sc_seed,1,1,8)
        self._L(fm,"Size ±%",1,2); self._E(fm,self.sc_spread,1,3,6)
        for k,ax in enumerate("xyz"):
            self._L(fm,f"{ax}0 / {ax}1",2+k,0)
            self._E(fm,self.sc_box[2*k],2+k,1,8); self._E(fm,self.sc_box[2*k+1],2+k,2,8)
        self._btn(fm,"whole device",self._def_sc_whole,bg=C["border"]).grid(row=2,column=3,columnspan=2,padx=2,sticky="w")
        bf=tk.Frame(fr3,bg=C["panel"]); bf.pack(fill=tk.X,padx=4,pady=2)
        self._btn(bf,"⁘ Scatter the form's defect",self._def_scatter,bg="#7a3fc4").pack(side=tk.LEFT,padx=2)
        gf=tk.Frame(fr3,bg=C["panel"]); gf.pack(fill=tk.X,padx=4,pady=2)
        tk.Label(gf,text="Group",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(4,2))
        self.sc_grp_cb=ttk.Combobox(gf,textvariable=self.sc_grp,values=[],state="readonly",width=28,font=self.F["ui"])
        self.sc_grp_cb.pack(side=tk.LEFT,padx=2)
        self._btn(gf,"✕ Remove group",self._def_rmgroup,bg="#a3243b").pack(side=tk.LEFT,padx=4)
        self._hint(fr3,"Scatters copies of the defect in the form at random centres inside the box (reproducible "
                       "with the seed; sizes varied by ±%). Density: per cm³ of the box, or per cm² of the plane "
                       "normal to the axis for lines through the whole device (e.g. a threading-dislocation "
                       "density). Scattered defects are ordinary rows afterwards, kept together as a group.")
        self._hint(fr3,f"Mesh: the graded mesh is made denser through the first {SD.REFINE_MAX} voids, particles, "
                       "lines and planes (the node count stays Nx × Ny × Nz; boxes get exact faces). Numbers of "
                       "states are exact on any mesh - deep levels in a shape smaller than the cells are spread "
                       "over the nearest nodes, and interface states carry the true surface area. A void or "
                       "particle smaller than the cells around it cannot replace whole cells without growing "
                       "several times its size, so it acts through its surface states only; the summary above "
                       "says so after a run - raise Nx/Ny/Nz (or add dense planes) to resolve it.")
        for e in self.df_ce+self.df_se+[self.df_dens_e,self.df_dit_e,self.df_q_e]+self.df_lv_w[2:]:
            e.configure(disabledbackground="#0b1422",disabledforeground="#46566b")
        self._def_cat_changed(); self._def_kind_changed()
        self._def_sc_unit()

    # ─── DEFECT CRUD ───────────────────────────────────────────
    def _def_cat_changed(self):
        names=[n for n,d in SD.LIB.items() if d["cat"]==self.df_cat.get()]
        self.df_kind_cb.configure(values=names)
        if names and self.df_kind.get() not in names:
            self.df_kind.set(names[0]); self._def_kind_changed()
    def _def_kind_changed(self,fill=True):
        p=SD.LIB.get(self.df_kind.get())
        if p is None: return
        self.df_desc.config(text=p["desc"])
        if fill:
            d=SD.Defect(self.df_kind.get())
            self._def_to_form(d,label=False,centre=False)
            self._def_fill_level(p)
        self._def_fields()
    def _def_fill_level(self,p):
        if p["levels"]:
            lv=p["levels"][0]
            self.df_lt.set(lv.type); self.df_lr.set(lv.ref); self.df_lE.set(f"{lv.E:g}")
            self.df_lsn.set(f"{lv.sn:g}"); self.df_lsp.set(f"{lv.sp:g}")
    def _def_fields(self):
        """Enable the form fields the chosen defect uses; units and size labels."""
        p=SD.LIB.get(self.df_kind.get())
        if p is None: return
        g=p["geom"]
        shapes={"inclusion":["sphere","box","cylinder"],"line":["cylinder"],"plane":["slab"]}.get(g,SD.SHAPES)
        self.df_shape_cb.configure(values=shapes)
        if self.df_shape.get() not in shapes: self.df_shape.set(shapes[0])
        sh=self.df_shape.get()
        def st(w,on,ro=False): w.configure(state=("readonly" if ro else "normal") if on else "disabled")
        labs=SD.SIZE_LABELS.get(sh,("","",""))
        for e,l in zip(self.df_se,labs): st(e,bool(l))
        self.df_sz_txt.config(text=" / ".join(l for l in labs if l) if any(labs) else "(whole host material)")
        for e in self.df_ce: st(e,sh!="everywhere")
        st(self.df_axis_cb,sh in ("cylinder","slab"),True)
        dens_on=g in ("volume","uniform","fluence","line") or (g=="plane" and bool(p["levels"]))
        st(self.df_dens_e,dens_on)
        self.df_dens_u.config(text=SD.DENSITY_UNIT.get(g,"") if dens_on else "")
        st(self.df_dit_e,(g=="inclusion" and not p["pin"]) or g=="interface" or (g=="plane" and p["dit"]>0))
        q_on=g in ("interface","charge"); st(self.df_q_e,q_on)
        self.df_q_u.config(text=("q/cm²" if g=="interface" else "q/cm³") if q_on else "")
        st(self.df_metal_cb,bool(p["pin"]),True)
        for w in self.df_lv_w: st(w,bool(p["levels"]),isinstance(w,ttk.Combobox))
        if p["levels"]:
            unit={"fluence":" cm⁻¹ (introduction rate)","interface":" cm⁻²"}.get(g,"")
            txt=["Levels:"]+[f"  {lv.text()}"+(f", {'η' if g=='fluence' else '×'} {lv.w:g}{unit}"
                                                if lv.w!=1. else "") for lv in p["levels"]]
            if len(p["levels"])>1: txt.append("  (the level fields edit the first one)")
            self.df_lvl_txt.config(text="\n".join(txt))
        else:
            self.df_lvl_txt.config(text="")
    def _def_to_form(self,d,label=True,centre=True):
        self.df_shape.set(d.shape); self.df_axis.set(d.axis)
        for v,x in (((self.df_x,d.x),(self.df_y,d.y),(self.df_z,d.z)) if centre else ())+((self.df_a,d.a),(self.df_b,d.b),(self.df_c,d.c)):
            v.set(f"{x:.6g}")
        self.df_dens.set(f"{d.density:g}"); self.df_dit.set(f"{d.dit:g}"); self.df_q.set(f"{d.charge:g}")
        if d.metal: self.df_metal.set(d.metal)
        self.df_host.set(d.host or "(any)")
        if label: self.df_lbl.set(d.label)
        else: self.df_lbl.set("")
        if not centre:                                # keep the form's centre if it lies in the device
            try:
                L3=(float(self.Lx.get()),float(self.Ly.get()),float(self.Lz.get()))
                if not all(0<=float(v.get())<=Lq for v,Lq in zip((self.df_x,self.df_y,self.df_z),L3)):
                    self._def_centre()
            except Exception: self._def_centre()
    def _def_centre(self):
        for v,L in ((self.df_x,self.Lx),(self.df_y,self.Ly),(self.df_z,self.Lz)):
            try: v.set(f"{0.5*float(L.get()):.6g}")
            except Exception: pass
    def _def_from_form(self):
        k=self.df_kind.get(); p=SD.LIB.get(k)
        if p is None: messagebox.showerror("Defect","Choose a defect from the library."); return None
        try:
            f=lambda v:float(str(v.get()).strip() or 0.)
            x,y,z,a,b,c=(f(v) for v in (self.df_x,self.df_y,self.df_z,self.df_a,self.df_b,self.df_c))
            dens,dit,q=f(self.df_dens),f(self.df_dit),f(self.df_q)
            E,sn,sp=f(self.df_lE),f(self.df_lsn),f(self.df_lsp)
        except Exception:
            messagebox.showerror("Defect","Centre, size, density, D_it, charge and level fields must be numbers."); return None
        sh=self.df_shape.get(); labs=SD.SIZE_LABELS.get(sh,("","",""))
        if labs[0] and a<=0: messagebox.showerror("Defect",f"The {labs[0]} must be larger than 0 nm."); return None
        if sh=="box" and min(b,c)<=0: messagebox.showerror("Defect","A box needs dx, dy and dz > 0 nm."); return None
        if min(b,c)<0 or dens<0 or dit<0:
            messagebox.showerror("Defect","Sizes, densities and D_it cannot be negative."); return None
        level=None
        if p["levels"]:
            lv=p["levels"][0]
            if sn<=0 or sp<=0: messagebox.showerror("Defect","Capture cross sections must be > 0."); return None
            if (self.df_lr.get()!="Ei" and not 0.<=E<=6.5) or (self.df_lr.get()=="Ei" and abs(E)>3.3):
                messagebox.showerror("Defect","The level must lie in the gap: 0…Eg below Ec / above Ev, "
                                              "or −Eg/2…+Eg/2 from midgap."); return None
            cur=(self.df_lt.get(),self.df_lr.get(),E,sn,sp)
            if ((cur[0],cur[1])!=(lv.type,lv.ref) or abs(E-lv.E)>1e-9 or abs(sn/lv.sn-1.)>1e-6
                    or abs(sp/lv.sp-1.)>1e-6):
                level=cur
        host=self.df_host.get(); host="" if host=="(any)" else host
        return SD.Defect(k,shape=sh,x=x,y=y,z=z,a=a,b=b,c=c,axis=self.df_axis.get(),density=dens,dit=dit,
                         charge=q,metal=self.df_metal.get() if p["pin"] else p["metal"],level=level,host=host,
                         label=self.df_lbl.get().strip())
    def _adddef(self):
        d=self._def_from_form()
        if d: self.defects.append(d); self._rfd(); self._select(self.dtree,len(self.defects)-1,False)
    def _editdef(self):
        i=self._sel_idx(self.dtree)
        if i is None: messagebox.showinfo("Defect","Select a defect in the table first."); return
        d=self._def_from_form()
        if d: d.group=self.defects[i].group; self.defects[i]=d; self._rfd()
    def _dupdef(self):
        i=self._sel_idx(self.dtree)
        if i is None: return
        d=SD.Defect.from_dict(self.defects[i].to_dict()); d.label+=" copy"; d.group=""
        self.defects.insert(i+1,d); self._rfd(); self._select(self.dtree,i+1,False)
    def _loaddef(self):
        if getattr(self,"_noload",False): return
        i=self._sel_idx(self.dtree)
        if i is None or i>=len(self.defects): return
        d=self.defects[i]; p=d.preset()
        self.df_cat.set(p["cat"])
        self.df_kind_cb.configure(values=[n for n,q in SD.LIB.items() if q["cat"]==p["cat"]])
        self.df_kind.set(d.kind); self._def_kind_changed(fill=False)
        self._def_to_form(d)
        lv=d.levels()
        if lv:
            self.df_lt.set(lv[0].type); self.df_lr.set(lv[0].ref); self.df_lE.set(f"{lv[0].E:g}")
            self.df_lsn.set(f"{lv[0].sn:g}"); self.df_lsp.set(f"{lv[0].sp:g}")
        self._def_fields()
    def _deldef(self):
        i=self._sel_idx(self.dtree)
        if i is not None: del self.defects[i]; self._rfd()
    def _defup(self):
        i=self._sel_idx(self.dtree)
        if i: self.defects[i-1],self.defects[i]=self.defects[i],self.defects[i-1]; self._rfd(); self._select(self.dtree,i-1,False)
    def _defdn(self):
        i=self._sel_idx(self.dtree)
        if i is not None and i<len(self.defects)-1:
            self.defects[i],self.defects[i+1]=self.defects[i+1],self.defects[i]; self._rfd(); self._select(self.dtree,i+1,False)
    def _clrdef(self):
        if self.defects and messagebox.askyesno("Defects","Delete all defects?"): self.defects.clear(); self._rfd()
    def _def_sc_unit(self):
        self.sc_unit.config(text="defects" if self.sc_mode.get()=="count" else "cm⁻³ (lines: cm⁻²)")
    def _def_sc_whole(self):
        for k,L in enumerate((self.Lx,self.Ly,self.Lz)):
            self.sc_box[2*k].set("0"); self.sc_box[2*k+1].set(f"{float(L.get()):g}")
    def _def_scatter(self):
        d=self._def_from_form()
        if d is None: return
        if d.shape=="everywhere":
            messagebox.showinfo("Scatter","This defect fills its host material already - choose a shape to scatter it."); return
        try:
            L3=(float(self.Lx.get()),float(self.Ly.get()),float(self.Lz.get()))
            bx=[float(v.get()) for v in self.sc_box]; val=float(self.sc_n.get())
            seed=int(float(self.sc_seed.get())); spread=float(self.sc_spread.get())/100.
        except Exception:
            messagebox.showerror("Scatter","Number/density, seed, spread and box must be numbers."); return
        if bx[1]<=bx[0] and bx[3]<=bx[2] and bx[5]<=bx[4]: bx=[0.,L3[0],0.,L3[1],0.,L3[2]]
        lo=np.clip([bx[0],bx[2],bx[4]],0.,L3); hi=np.clip([bx[1],bx[3],bx[5]],0.,L3)
        if np.any(hi<=lo): messagebox.showerror("Scatter","The box must have x1 > x0, y1 > y0, z1 > z0 inside the device."); return
        if not 0.<=spread<1.: messagebox.showerror("Scatter","Size spread must be 0…99 %."); return
        if self.sc_mode.get()=="count": n=int(round(val)); how=f"{n}"
        else:
            nexp=SD.count_for_density(d,val,L3,(lo,hi)); n=int(round(nexp))
            how=f"{val:g} {'cm⁻²' if (d.shape=='cylinder' and d.b<=0) else 'cm⁻³'} → {nexp:.3g}"
            if n<1:
                messagebox.showinfo("Scatter",f"That density gives {nexp:.2g} defects in the box - fewer than one. "
                                               f"Enlarge the box or raise the density."); return
        if n<1: return
        if n>5000: messagebox.showerror("Scatter",f"{n} defects is too many (at most 5000)."); return
        if n>500 and not messagebox.askyesno("Scatter",f"Scatter {n} defects? Each one is mapped onto the mesh "
                                                      f"separately, which takes a few seconds per thousand."): return
        new=SD.scatter(d,n,(lo,hi),seed,spread)
        self.defects.extend(new); self._rfd()
        self.sc_grp.set(new[0].group)
        self._log(f"Scattered {how} × {d.kind} (seed {seed}) in x {lo[0]:g}-{hi[0]:g}, y {lo[1]:g}-{hi[1]:g}, "
                  f"z {lo[2]:g}-{hi[2]:g} nm")
    def _def_rmgroup(self):
        g=self.sc_grp.get()
        if not g: return
        n=sum(1 for d in self.defects if d.group==g)
        self.defects=[d for d in self.defects if d.group!=g]; self._rfd()
        self._log(f"Removed the {n} defects of group '{g}'")
    def _rfd(self,redraw=True):
        tv=self.dtree; sel=tv.selection(); idx=tv.index(sel[0]) if sel else None
        tv.delete(*tv.get_children())
        for k,d in enumerate(self.defects):
            cen="-" if d.shape=="everywhere" else f"{d.x:.4g}, {d.y:.4g}, {d.z:.4g}"
            tv.insert("","end",values=(k,d.label,d.kind,d.describe(),cen,d.strength_text(),d.host or "any"))
        if idx is not None and idx<len(self.defects): self._select(tv,idx,load=False)
        groups=[]
        for d in self.defects:
            if d.group and d.group not in groups: groups.append(d.group)
        self.sc_grp_cb.configure(values=groups)
        if self.sc_grp.get() not in groups: self.sc_grp.set(groups[-1] if groups else "")
        self._def_summary()
        if redraw: self._geometry_changed()
    def _def_summary(self,notes=None):
        if not self.defects: txt="no defects"
        else:
            n=len(self.defects); ng=len({d.group for d in self.defects if d.group})
            nm=sum(1 for d in self.defects if d.preset()["geom"] in ("inclusion","line","plane") and d.shape!="everywhere")
            txt=f"{n} defect{'s' if n>1 else ''}"+(f" ({ng} scattered group{'s' if ng>1 else ''})" if ng else "")
            if nm>SD.REFINE_MAX: txt+=f" · {nm-SD.REFINE_MAX} voids/particles/lines/planes beyond the first {SD.REFINE_MAX} get no mesh refinement"
        if notes: txt+="\n⚠ "+"\n⚠ ".join(notes)
        self.def_sum.config(text=txt,fg=C["yellow"] if notes else C["green"])

    # ─── DEFECT DRAWING ────────────────────────────────────────
    _DEF_DRAW_MAX=1500
    _DEF_CAT_COL={"Point defects & contamination":"#fb923c","Radiation damage":"#e879f9","Interfaces":"#22d3ee",
                  "Extended defects":"#a3e635","SiC · GaN · GaAs":"#f472b6"}
    @classmethod
    def _def_color(cls,d):
        p=d.preset()
        if p["geom"]=="inclusion":
            return "#facc15" if p["pin"] else ("#93c5fd" if p["material"]=="SiO2" else "#f1f5f9")
        return cls._DEF_CAT_COL.get(p["cat"],"#fb923c")
    def _draw_defects_2d(self,ax,defects,H,Vv,N,pos,L3,overlay=False):
        """Defects in a section plane: those the plane cuts filled (overlay on a
        solution map: outlined), the others as dotted outlines of their
        projection (structure view only).  Returns the number cut."""
        from matplotlib.patches import Circle,Rectangle
        if not defects: return 0
        span=max(L3["xyz".index(H)],L3["xyz".index(Vv)]); thin=0.004*span; tol=1e-9*max(L3); ncut=0
        for d in defects[:self._DEF_DRAW_MAX]:
            sec=defect_section(d,H,Vv,N,pos,L3,tol)
            if sec is None or (overlay and not sec[0]): continue
            cut,kind,g=sec; col=self._def_color(d); ncut+=bool(cut)
            if cut and not overlay: kw=dict(facecolor=col,edgecolor=col,alpha=0.6,lw=1.0,zorder=5)
            elif cut: kw=dict(facecolor="none",edgecolor=col,alpha=0.95,lw=1.1,zorder=6)
            else: kw=dict(facecolor="none",edgecolor=col,alpha=0.5,lw=0.9,ls=":",zorder=4)
            if kind=="circle":
                ax.add_patch(Circle((g[0],g[1]),max(g[2],0.5*thin),**kw))
            else:
                h0,v0,h1,v1=g
                if h1-h0<thin: m=0.5*(h0+h1); h0,h1=m-0.5*thin,m+0.5*thin
                if v1-v0<thin: m=0.5*(v0+v1); v0,v1=m-0.5*thin,m+0.5*thin
                ax.add_patch(Rectangle((h0,v0),h1-h0,v1-v0,**kw))
        return ncut
    def _draw_defects_3d(self,a3,defects,L3):
        if not defects: return
        xs,zs,ys,ss,cs=[],[],[],[],[]; Lm=max(L3)
        for d in defects[:self._DEF_DRAW_MAX]:
            if d.shape=="everywhere": continue
            col=self._def_color(d); lo,hi=d.extent(L3)
            if d.shape in ("sphere","box"):
                xs.append(d.x); zs.append(d.z); ys.append(-d.y); cs.append(col)
                ss.append(float(np.clip((float(np.max(hi-lo))/Lm*260.)**2,6.,500.)))
            elif d.shape=="cylinder":
                q="xyz".index(d.axis); a=np.array([d.x,d.y,d.z]); b=a.copy(); a[q]=lo[q]; b[q]=hi[q]
                a3.plot([a[0],b[0]],[a[2],b[2]],[-a[1],-b[1]],color=col,lw=1.6,alpha=0.9)
            else:
                q="xyz".index(d.axis); o=[k for k in range(3) if k!=q]; pts=[]
                for u,v in ((0,0),(1,0),(1,1),(0,1)):
                    pp=np.array([d.x,d.y,d.z],float); pp[o[0]]=(lo,hi)[u][o[0]]; pp[o[1]]=(lo,hi)[v][o[1]]
                    pts.append((pp[0],pp[2],-pp[1]))
                a3.add_collection3d(Poly3DCollection([pts],alpha=0.22,facecolor=col,edgecolor=col,linewidth=0.5))
        if xs: a3.scatter(xs,zs,ys,s=ss,c=cs,depthshade=False,edgecolors="#111111",linewidths=0.3,alpha=0.95)
    def _defect_legend(self,defects):
        from matplotlib.lines import Line2D
        if not defects: return []
        by={}
        for d in defects: by.setdefault(d.kind,[]).append(d)
        out=[Line2D([0],[0],color="none",label="Defects (dotted: off this plane)")]
        for k,ds in by.items():
            d=ds[0]; col=self._def_color(d)
            where=(" everywhere"+(f" in {d.host}" if d.host else "")) if d.shape=="everywhere" else ""
            lab=f"{k}{f' ×{len(ds)}' if len(ds)>1 else ''}: {d.strength_text()}{where}"
            out.append(Line2D([0],[0],marker="o",ls="none",markersize=7,markerfacecolor=col,markeredgecolor=col,
                              label=textwrap.fill(lab,34,subsequent_indent="  ")))
        return out
    def _defect_at(self,P,H,Vv,N,pos,L3):
        """The defect drawn under the point P (dict of x, y, z) in the section, or None."""
        span=max(L3["xyz".index(H)],L3["xyz".index(Vv)]); thin=0.006*span
        for d in reversed(self.defects[:self._DEF_DRAW_MAX]):
            sec=defect_section(d,H,Vv,N,pos,L3)
            if sec is None or not sec[0]: continue
            _,kind,g=sec
            if kind=="circle":
                if (P[H]-g[0])**2+(P[Vv]-g[1])**2<=max(g[2],thin)**2: return d
            elif min(g[0],g[2])-thin<=P[H]<=max(g[0],g[2])+thin and min(g[1],g[3])-thin<=P[Vv]<=max(g[1],g[3])+thin:
                return d
        return None

    def _page_solver(self,p):
        fr=self._sec(p,"Physics")
        fm=tk.Frame(fr,bg=C["panel"]); fm.pack(fill=tk.X,padx=4,pady=4)
        self._L(fm,"Temperature (K)",0,0); self._E(fm,self.T,0,1,8)
        self._L(fm,"Background SRH τ0 (s)",10,0); self._E(fm,self.tau0,10,1,8)
        self._L(fm,"empty = material",10,2)
        for k,(txt,var) in enumerate([("Band-gap narrowing (Klaassen)",self.m_bgn),
                                      ("SRH lifetime τ(N)",self.m_tau),
                                      ("Velocity saturation",self.m_vsat),
                                      ("Quantum confinement at oxide interfaces (MLDA)",self.m_qc),
                                      ("Surface mobility (Lombardi, Si)",self.m_lb),
                                      ("Gate tunnelling (direct + Fowler-Nordheim)",self.m_tun),
                                      ("Metal contacts: thermionic emission + image-force lowering",self.m_te),
                                      ("Tunnelling through Schottky barriers (field emission)",self.m_sbt),
                                      ("Impact ionisation (avalanche breakdown)",self.m_ii)]):
            ttk.Checkbutton(fm,text=txt,variable=var).grid(row=1+k,column=0,columnspan=3,sticky="w",padx=4,pady=1)
        self._hint(fr,"MLDA, Lombardi and gate tunnelling matter for thin-oxide MOS devices: MLDA pushes the "
                      "inversion layer ~1 nm off the oxide (V_th up, gate capacitance down), Lombardi lowers "
                      "the channel mobility with the normal field, and tunnelling lets current flow through "
                      "the gate stack (see each gate's stack on the Contacts page). The metal-contact models act "
                      "on Schottky contacts and on contacts with a metal: carriers cross the interface by "
                      "thermionic emission over the image-force-lowered barrier and by tunnelling through it "
                      "(Tsu-Esaki with WKB), so a metal on heavily doped silicon is ohmic and on lightly doped "
                      "silicon rectifies. Off: a Schottky contact is the classic fixed-density barrier and an "
                      "Ohmic contact with a metal stays ideal. Templates set these switches.")
        self._hint(fr,"Impact ionisation: carriers accelerated by a high field knock out electron-hole pairs "
                      "(Chynoweth law with the van Overstraeten-de Man coefficients for Si; GaAs, 4H-SiC and GaN "
                      "from the literature, see the Help). It turns a reverse-biased junction's current into an "
                      "avalanche: the I-V tab then reports the breakdown voltage where the current runs away. "
                      "Near breakdown the multiplication grows without bound and the solver stops converging; "
                      "the last converged points mark BV. Off: fields can be read, breakdown cannot.")
        self._hint(fr,"Background SRH τ0: the lifetime of the semiconductors' built-in midgap recombination "
                      "(empty = material defaults: Si 10 µs, Ge 1 ms, 4H-SiC 1 µs for electrons and 0.1 µs for "
                      "holes, GaN 10 ns, GaAs 1 ns; τ(N) still "
                      "divides it by 1 + N/Nref). 1e-3 = float-zone-quality silicon, 'none' = no background SRH at "
                      "all, so only the deep levels of the Defects page recombine - the clean way to study what "
                      "a defect does.")
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
        if FLAT_BUTTONS:       # macOS draws tk.Checkbutton natively: use the themed one
            ttk.Checkbutton(parent,text=text,variable=var,command=cmd,
                            style=self._check_style(fg or C["text"])).pack(side=tk.LEFT,padx=4)
            return
        tk.Checkbutton(parent,text=text,variable=var,command=cmd,bg=C["panel"],fg=fg or C["text"],
                       selectcolor="#07111e",activebackground=C["panel"],activeforeground=fg or C["text"],
                       font=self.F["ui"],highlightthickness=0).pack(side=tk.LEFT,padx=4)
    def _check_style(self,fg):
        st=f"C{fg.lstrip('#')}.TCheckbutton"
        done=self.__dict__.setdefault("_chk_styles",set())
        if st not in done:
            s=ttk.Style(self)
            s.configure(st,foreground=fg,background=C["panel"],font=self.F["ui"])
            s.map(st,foreground=[("disabled",C["sub"]),("active",fg)],background=[("active",C["panel"])])
            done.add(st)
        return st

    # ─── PANEL VISIBILITY (Panels ▾ menu of every plot tab) ───
    _PANELS={"td":[("3d","3-D view"),("sec","Cross-section"),("leg","Legend")],
             "ts":[("3d","3-D slices"),("xy","XY plane"),("yz","YZ plane"),("xz","XZ plane")],
             "tb":[("band","Band diagram"),("ribbon","Band ribbon (3-D)")],
             "tc":[("line","Line cut: current density"),("map","In-plane map"),
                   ("prof","Current under a contact"),("rec","Line cut: recombination"),("3d","3-D current field")],
             "tf":[("E","Electric field"),("mu","Mobility")],
             "ti":[("lin","Linear I-V"),("log","Semi-log I-V (all terminals)"),("bjt","BJT: gain and f_T vs I_C")]}
    def _panels_menu(self,parent,tab,cmd,side=tk.RIGHT):
        """A 'Panels ▾' menu button: show / hide each panel of a plot tab
        (the panels left take the room of the hidden ones)."""
        mb=ttk.Menubutton(parent,text="Panels ▾",style="TMenubutton")
        m=tk.Menu(mb,tearoff=0,bg=C["panel"],fg=C["text"],activebackground=C["accent"],
                  activeforeground="white",selectcolor=C["green"],font=self.F["ui"])
        d=self.pv.setdefault(tab,{})
        for key,label in self._PANELS[tab]:
            v=d.setdefault(key,tk.BooleanVar(value=True))
            m.add_checkbutton(label=label,variable=v,command=cmd)
        m.add_separator()
        m.add_command(label="Show all",command=lambda:(self._panels_all(tab,True),cmd()))
        mb["menu"]=m; mb.pack(side=side,padx=4)
        Tip(mb,"Show or hide the panels of this tab - the panels left share the room.",self.F["small"])
        return mb
    def _panels_all(self,tab,on):
        for v in self.pv.get(tab,{}).values(): v.set(on)
    def _vis(self,tab):
        """Keys of the visible panels of a tab, in display order."""
        d=self.pv.get(tab,{})
        return [k for k,_ in self._PANELS[tab] if (d[k].get() if k in d else True)]
    def _no_panels(self,fig,cv):
        self._placeholder(fig,cv,"All panels of this tab are hidden - use Panels ▾ to show them")

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
        lb=self._btn(cf,"Layers…",self._layers_win,bg=C["border"])
        lb.pack(side=tk.LEFT,padx=4)
        Tip(lb,"Show or hide regions, contacts and defects in the 3-D view (and optionally the\n"
               "cross-section), make contacts see-through, show the insulators in 3-D.",self.F["small"])
        self._panels_menu(cf,"td",self._ddev)
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
              "log10|J gate|","Λn quantum","trap density","trapped charge","log10 G avalanche"]
        cb=ttk.Combobox(cf,textvariable=self.sl_qty,values=qtys,state="readonly",width=13,font=self.F["ui"])
        cb.pack(side=tk.LEFT,padx=2); cb.bind("<<ComboboxSelected>>",lambda _:self._upd_slice())
        ttk.Separator(cf,orient="vertical").pack(side=tk.LEFT,fill=tk.Y,padx=6,pady=2)
        self._check(cf,"Regions",self.show_regions,self._redraw_solution,"#f0abfc")
        self._check(cf,"Labels",self.show_labels,self._redraw_solution,"#f0abfc")
        self._panels_menu(cf,"ts",self._upd_slice,side=tk.LEFT)
        self.sl_hover=tk.Label(cf,text="",bg=C["panel"],fg=C["text"],font=self.F["small"],anchor="w")
        self.sl_hover.pack(side=tk.LEFT,fill=tk.X,expand=True,padx=(10,4))
        sf=self._strip(self.ts)
        tk.Label(sf,text="Position",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=4)
        self.sl_scale=ttk.Scale(sf,from_=0,to=30,orient=tk.HORIZONTAL,variable=self.sl_pos,
                                command=lambda _:self._debounce('_sl_after_id',self._upd_slice))
        self.sl_scale.pack(side=tk.LEFT,fill=tk.X,expand=True,padx=6)
        self.sl_lbl=tk.Label(sf,text="",bg=C["panel"],fg=C["yellow"],font=self.F["mono"],width=22,anchor="w")
        self.sl_lbl.pack(side=tk.LEFT,padx=4)
        self.fsl,self.cvsl=self._mfig(self.ts,fs=(15,7))
        self.cvsl.mpl_connect("motion_notify_event",lambda e:self._map_hover(e,self.sl_hover))
        self.cvsl.mpl_connect("draw_event",self._labels_on_draw)

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
            self._check(cf,"Labels",self.show_labels,self._redraw_solution,"#f0abfc")
            self._panels_menu(cf,"tb",self._upd_band_plot,side=tk.LEFT)
        self.bd_lbl1,self.bd_lbl2,self.bd_scale1,self.bd_scale2=self._build_cut_controls(
            self.tb,self.bd_axis,self.bd_p1,self.bd_p2,self._upd_band_axis,
            lambda _:self._debounce('_bd_after_id',self._upd_band_plot),extras)
        self.fb,self.cvb=self._mfig(self.tb)

    def _arrow_n(self):
        """Arrow counts (2-D: across the map, 3-D: along x), clipped to their ranges."""
        out=[]
        for v,(lo,hi),d in ((self.cur_n2d,ARROWS_2D,22),(self.cur_n3d,ARROWS_3D,8)):
            try: out.append(max(lo,min(hi,int(v.get()))))
            except Exception: out.append(d)
        return tuple(out)
    def _arrows_changed(self,_=None):
        n2,n3=self._arrow_n()
        if (self.cur_n2d.get() if self._intok(self.cur_n2d) else None)!=n2: self.cur_n2d.set(n2)
        if (self.cur_n3d.get() if self._intok(self.cur_n3d) else None)!=n3: self.cur_n3d.set(n3)
        self._debounce('_cur_after_id',self._upd_curr_plot,250)
    @staticmethod
    def _intok(v):
        try: v.get(); return True
        except Exception: return False
    def _arrow_spin(self,cf,var,rng,step,tip):
        sp=ttk.Spinbox(cf,from_=rng[0],to=rng[1],increment=step,width=4,textvariable=var,
                       command=self._arrows_changed,font=self.F["ui"])
        sp.pack(side=tk.LEFT,padx=(0,6))
        for ev in ("<Return>","<KP_Enter>","<FocusOut>"): sp.bind(ev,self._arrows_changed,add="+")
        Tip(sp,tip,self.F["small"])
        return sp

    def _init_curr_tab(self):
        def extras(cf):
            self._check(cf,"2-D arrows",self.cur_show_arrows2d,self._upd_curr_plot,C["green"])
            self._arrow_spin(cf,self.cur_n2d,ARROWS_2D,2,
                             "Number of arrows across the in-plane map (the rows follow the map's shape).\n"
                             "More arrows show finer detail of the current paths.")
            self._check(cf,"3-D arrows",self.cur_show_arrows3d,self._upd_curr_plot,C["yellow"])
            self._arrow_spin(cf,self.cur_n3d,ARROWS_3D,1,
                             "Number of arrows along x in the 3-D field; depth y gets 3/4 and z 5/8 of it\n"
                             "(never more than the mesh lines along an axis, so a 2-D cross-section keeps\n"
                             f"its few z planes; at most {ARROWS_3D_MAX} arrows in all).")
            self._check(cf,"Regions",self.show_regions,self._redraw_solution,"#f0abfc")
            self._check(cf,"Labels",self.show_labels,self._redraw_solution,"#f0abfc")
            self.cur_hover=tk.Label(cf,text="",bg=C["panel"],fg=C["text"],font=self.F["small"],anchor="w")
            self.cur_hover.pack(side=tk.LEFT,fill=tk.X,expand=True,padx=(10,4))
            cf2=self._strip(self.tc)                  # second row: contact profile, structure, panels
            tk.Label(cf2,text="Contact",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(4,2))
            self.cb_prof=ttk.Combobox(cf2,textvariable=self.cur_prof,values=["(none)"],state="readonly",width=22,
                                      font=self.F["ui"])
            self.cb_prof.pack(side=tk.LEFT,padx=2); self.cb_prof.bind("<<ComboboxSelected>>",self._prof_pick)
            Tip(self.cb_prof,"Current density through a contact, node by node (the 'Current under a contact'\n"
                             "panel): how evenly it carries its current. A wide emitter crowds its current\n"
                             "to the edges nearest the base contact - the base current's lateral voltage\n"
                             "drop turns the middle off. Face contacts only.",self.F["small"])
            ttk.Separator(cf2,orient="vertical").pack(side=tk.LEFT,fill=tk.Y,padx=6,pady=2)
            tk.Label(cf2,text="Map",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(2,2))
            cb=ttk.Combobox(cf2,textvariable=self.cur_bg,values=["|J|","structure","|J| + structure"],state="readonly",
                            width=16,font=self.F["ui"])
            cb.pack(side=tk.LEFT,padx=2); cb.bind("<<ComboboxSelected>>",lambda e:self._upd_curr_plot())
            Tip(cb,"What lies under the arrows of the in-plane map: the current density |J|, the device\n"
                   "structure (region colours, see-through), or both.",self.F["small"])
            self._check(cf2,"Structure in 3-D",self.cur_st3,self._upd_curr_plot,"#f0abfc")
            tk.Label(cf2,text="Opacity",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(6,2))
            sc=ttk.Scale(cf2,from_=0.,to=1.,orient=tk.HORIZONTAL,variable=self.cur_sa,length=120,
                         command=lambda *_:self._debounce("_cur_sa_id",self._upd_curr_plot,200))
            sc.pack(side=tk.LEFT,padx=2)
            self._cur_sa_lbl=tk.Label(cf2,text=f"{self.cur_sa.get():.2f}",width=4,bg=C["panel"],fg=C["text"],font=self.F["ui"])
            self._cur_sa_lbl.pack(side=tk.LEFT)
            def _sa_txt(*_):
                try: self._cur_sa_lbl.config(text=f"{float(self.cur_sa.get()):.2f}")
                except (tk.TclError,ValueError): pass
            self.cur_sa.trace_add("write",_sa_txt)      # also when set from a script or a loaded file
            Tip(sc,"Opacity of the structure drawn with the current arrows (map and 3-D). Parts\n"
                   "hidden in Layers… stay hidden; a part given its own opacity there keeps it.",self.F["small"])
            b=self._btn(cf2,"Layers…",self._layers_win,bg=C["border"]); b.pack(side=tk.LEFT,padx=6)
            self._panels_menu(cf2,"tc",self._upd_curr_plot,side=tk.LEFT)
        self.cur_lbl1,self.cur_lbl2,self.cur_scale1,self.cur_scale2=self._build_cut_controls(
            self.tc,self.cur_axis,self.cur_p1,self.cur_p2,self._upd_curr_axis,
            lambda _:self._debounce('_cur_after_id',self._upd_curr_plot),extras)
        self.fc,self.cvc=self._mfig(self.tc)
        self.cvc.mpl_connect("motion_notify_event",lambda e:self._map_hover(e,self.cur_hover))
        self.cvc.mpl_connect("draw_event",self._labels_on_draw)

    def _prof_choices(self):
        cb=getattr(self,"cb_prof",None)
        if cb is None: return
        vals=["(none)"]+list(self._gmap.keys()); cb.configure(values=vals)
        if self.cur_prof.get() not in vals: self.cur_prof.set("(none)")
    def _prof_pick(self,_=None):
        if self.cur_prof.get()!="(none)": self.pv["tc"]["prof"].set(True)
        self._upd_curr_plot()

    def _init_fi_tab(self):
        def extras(cf):
            self._check(cf,"E-field",self.fi_show_field,self._upd_fi_plot,C["accent"])
            self._check(cf,"Mobility",self.fi_show_mob,self._upd_fi_plot,C["green"])
            self._check(cf,"Regions",self.show_regions,self._redraw_solution,"#f0abfc")
            self._check(cf,"Labels",self.show_labels,self._redraw_solution,"#f0abfc")
            self._panels_menu(cf,"tf",self._upd_fi_plot,side=tk.LEFT)
        self.fi_lbl1,self.fi_lbl2,self.fi_scale1,self.fi_scale2=self._build_cut_controls(
            self.tf,self.fi_axis,self.fi_p1,self.fi_p2,self._upd_fi_axis,
            lambda _:self._debounce('_fi_after_id',self._upd_fi_plot),extras)
        self.ff,self.cvf=self._mfig(self.tf)

    def _init_iv_tab(self):
        cf=self._strip(self.ti)
        b=self._btn(cf,"Pin sweep",self._pin_sweep,bg=C["border"]); b.pack(side=tk.LEFT,padx=4)
        Tip(b,"Keep the last sweep as a reference curve: later sweeps are drawn over it, e.g.\n"
              "before and after a design change. Pins… lists, hides, renames and removes them.",self.F["small"])
        self._btn(cf,"Pins…",self._pins_win,bg=C["border"]).pack(side=tk.LEFT,padx=4)
        b=self._btn(cf,"A/B compare…",self._ab_win,bg=C["border"]); b.pack(side=tk.LEFT,padx=4)
        Tip(b,"Capture two designs - two templates, two variants of one, or one device before\n"
              "and after an edit - and run the same sweep on both: their curves come out side\n"
              "by side, each pinned, with the parameters of both.",self.F["small"])
        b=self._btn(cf,"Table…",self._ptab_win,bg=C["border"]); b.pack(side=tk.LEFT,padx=4)
        Tip(b,"The parameters of the last sweep and of every pinned sweep side by side - V_th, SS,\n"
              "I_on/I_off, n, V_A, BV, β, f_T, V_M … - with their differences from a reference\n"
              "column; copy it into a spreadsheet or export it as CSV.",self.F["small"])
        self.iv_pin_lbl=tk.Label(cf,text="",bg=C["panel"],fg=C["sub"],font=self.F["small"],anchor="w")
        self.iv_pin_lbl.pack(side=tk.LEFT,padx=8)
        self._panels_menu(cf,"ti",self._replot_iv)
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
        tk.Label(mid,text="Extracted parameters (last sweep and pinned sweeps)",bg=C["panel"],fg=C["accent"],
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
        self._ptab_fill()

    # ─── TEMPLATE LOADING ──────────────────────────────────────
    def _load(self,name):
        if self.busy:
            messagebox.showinfo("Busy","A solve is running - stop it first (Esc)."); return
        try: t=TEMPLATES[name]()
        except Exception:
            messagebox.showerror("Template error",traceback.format_exc()); return
        self._apply_template(name,t)
    def _apply_template(self,name,t,base=None):
        """Put a template dict into the editor; base = the template a variant
        belongs to (its documentation and menu category)."""
        base=base or name
        self.tname=name
        self.tmeta={k:copy.deepcopy(t.get(k)) for k in ("nodes","clear","refine","view","ivnet","ivfloat","ivout","wl","extrude")}
        self.mesh_type.set(t.get("mesh","rect")); self.prism_ax.set("auto")
        self.Lx.set(t["Lx"]); self.Ly.set(t["Ly"]); self.Lz.set(t["Lz"])
        self.Nx.set(t["Nx"]); self.Ny.set(t["Ny"]); self.Nz.set(t["Nz"])
        self.use_graded.set(True)
        self.regs=t["regs"]; self.cons=t["cons"]; self.sheets=t.get("sheets",[])
        self.defects=[SD.Defect.from_dict(d.to_dict()) for d in (t.get("defects") or [])]
        self.amr_nodes={}; self._amr_label()
        self.dense_x=list(t.get("dense_x",[])); self.dense_y=list(t.get("dense_y",[])); self.dense_z=list(t.get("dense_z",[]))
        self.ivS.set(t["ivS"]); self.ivE.set(t["ivE"]); self.ivN.set(t["ivN"])
        self.bv_thr.set(t.get("bv_thr",0.))
        self.tau0.set("" if t.get("tau0") is None else ("none" if t["tau0"]==0 else f"{t['tau0']:g}"))
        mdl=int(t.get("models",M_DEFAULT))
        for v,b in ((self.m_bgn,M_BGN),(self.m_tau,M_TAU),(self.m_vsat,M_VSAT),(self.m_qc,M_QC),
                    (self.m_lb,M_LB),(self.m_tun,M_TUN),(self.m_te,M_TE),(self.m_sbt,M_SBT),
                    (self.m_ii,M_II)): v.set(bool(mdl&b))
        groups=con_groups(self.cons)
        tgt=t.get("ivnet")
        if not tgt:
            ivc=t.get("ivc",0); tgt=next((k for k,idx in groups if ivc in idx),groups[0][0] if groups else "")
        self._refresh_groups(tgt=tgt,flt=t.get("ivfloat") or "(none)",out=t.get("ivout") or "auto")
        self._rfr(redraw=False); self._rfc(redraw=False); self._rfd(redraw=False); self._def_centre(); self._def_sc_whole()
        pk=t.get("prof")                  # Current tab: show the current under this contact
        if pk and group_of(self.cons,pk):
            disp=group_label(self.cons,pk,group_of(self.cons,pk))
            if disp in self._gmap: self.cur_prof.set(disp); self.pv["tc"]["prof"].set(True)
        cat=next((c for c,ns in template_groups() if base in ns),"")
        self.t_name.config(text=name); self.t_cat.config(text=cat)
        self.t_doc.configure(state="normal"); self.t_doc.delete("1.0","end")
        self.t_doc.insert("end",template_doc(base)); self.t_doc.configure(state="disabled")
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
                    defects=[SD.Defect.from_dict(d.to_dict()) for d in self.defects] if deep else self.defects,
                    dense_x=list(self.dense_x),dense_y=list(self.dense_y),dense_z=list(self.dense_z),
                    nodes=cp(self.tmeta.get("nodes")),clear=cp(self.tmeta.get("clear")),
                    refine=cp(self.tmeta.get("refine")),
                    mesh=self.mesh_type.get(),prism=None if self.prism_ax.get()=="auto" else self.prism_ax.get(),
                    extrude=self.tmeta.get("extrude"),tau0=self._tau0(strict=False),
                    amr=copy.deepcopy(self.amr_nodes) if self.amr_nodes else None)
    def _amr_clear(self):
        if self.busy: return
        had=bool(self.amr_nodes); self.amr_nodes={}; self._amr_label()
        if had: self._log("Adaptive mesh refinement cleared"); self._mesh_update()
    def _amr_label(self,info=None):
        n={k:len(v) for k,v in (self.amr_nodes or {}).items() if v}
        if not n: self.amr_lbl.config(text="no refinement"); return
        txt="refined: +"+", ".join(f"{len(self.amr_nodes[k])} {k}" for k in "xyz" if self.amr_nodes.get(k))+" lines"
        if info: txt+="  ("+"; ".join(f"pass {i+1}: +{a} → {b:,} nodes{'' if c else ', NOT converged'}"
                                     for i,(a,b,c) in enumerate(info))+")"
        self.amr_lbl.config(text=txt)
    def _tau0(self,strict=True):
        """Background SRH lifetime from the Solver page: None = material
        defaults, 0 = none.  strict: a bad entry raises ValueError."""
        v=self.tau0.get().strip().lower()
        if v in ("","default","material"): return None
        if v in ("none","off","inf","∞","0"): return 0.
        try:
            x=float(v)
            if x<0: raise ValueError
            return x
        except ValueError:
            if strict: raise ValueError("Background SRH τ0 (Solver page) must be a lifetime in s, 'none' or empty.")
            return None
    def _mesh(self):
        """Node coordinates (nm) of the present device, cached."""
        try: t=self._device_dict(deep=False)
        except Exception: return None
        key=repr((t["Lx"],t["Ly"],t["Lz"],t["Nx"],t["Ny"],t["Nz"],self.use_graded.get(),
                  [(r.mat,r.x0,r.y0,r.z0,r.x1,r.y1,r.z1,r.dtype,r.doping,_shape_key(r)) for r in t["regs"]],
                  [(c.face,c.i0p,c.i1p,c.j0p,c.j1p,c.k0p,c.k1p,c.bc) for c in t["cons"]],
                  [(s.axis,s.pos) for s in t["sheets"]],t["dense_x"],t["dense_y"],t["dense_z"],
                  repr(t["nodes"]),repr(t.get("clear")),repr(t["refine"]),
                  [(d.kind,d.shape,d.x,d.y,d.z,d.a,d.b,d.c,d.axis) for d in t.get("defects") or []
                   if d.preset()["geom"] in ("inclusion","line","plane")],repr(t.get("amr"))))
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
        self._lay_prune()
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
        L3=None
        for k,c in enumerate(self.cons):
            gate=c.bc==2 and c.face!=6
            tx=(f"{c.eot():.3g}"+(f" ({Con.stack_str(c.stack)})" if c.stack else "")) if gate else "-"
            pm=f"{c.pm:g}" if c.bc else "-"
            if getattr(c,"metal",""):
                try:
                    if L3 is None: L3=(float(self.Lx.get()),float(self.Ly.get()),float(self.Lz.get()))
                    bc_c,pm_c,_=contact_core(c,self.regs,L3)
                    pm=f"{pm_c:.3g}"
                except Exception: pass
            self.ctree.insert("","end",values=(k,c.label+(" ◯" if getattr(c,"shape",None) is not None else ""),
                c.net or "",Con.FNAME[c.face],Con.BCNAME[c.bc],getattr(c,"metal","") or "-",
                f"{c.V:g}",pm,tx,
                f"{c.nox:.1e}" if gate and c.nox else "-",self._con_range(c)))
        if idx is not None and idx<len(self.cons):
            self._select(self.ctree,idx,load=False)
        self._refresh_groups(); self._lay_prune(); self._prof_choices()
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
        metal=self.cb_metal.get(); metal=metal if metal in METALS else ""
        c=Con(face,i0/100.,i1/100.,j0/100.,j1/100.,bc,V,pm,lbl,tox,nox,k0/100.,k1/100.,
              net=self.cb_net.get().strip(),stack=stack,metal=metal)
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
        self.cb_stack.set(Con.stack_str(c.stack)); self.cb_metal.set(getattr(c,"metal","") or METAL_NONE)
    def _upd_metal(self):
        """Barrier of the metal chosen in the contact form, on whatever lies
        under the contact; the φm entry is locked while a metal sets it."""
        if not hasattr(self,"metal_lbl"): return
        metal=self.cb_metal.get()
        if metal not in METALS:
            self.metal_lbl.config(text=""); self.cb_pm_e.config(state="normal"); return
        FM={n:k for k,n in enumerate(Con.FNAME)}; BM={n:k for k,n in enumerate(Con.BCNAME)}
        try:
            c=Con(FM[self.cb_face.get()],float(self.cb_i0.get())/100.,float(self.cb_i1.get())/100.,
                  float(self.cb_j0.get())/100.,float(self.cb_j1.get())/100.,BM[self.cb_bc.get()],
                  k0p=float(self.cb_k0.get())/100.,k1p=float(self.cb_k1.get())/100.,metal=metal)
            L3=(float(self.Lx.get()),float(self.Ly.get()),float(self.Lz.get()))
            bc_c,pm_c,info=contact_core(c,self.regs,L3,self._models())
        except Exception:
            return
        self.cb_pm.set(round(pm_c,4)); self.cb_pm_e.config(state="disabled")
        self.metal_lbl.config(text=info+(f"  →  φm used {pm_c:.3f} eV" if bc_c==1 else ""))
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
            raise ValueError("Domain, mesh and solver fields must be numbers (background lifetime: "
                             "a time in s, 'none', or empty for the material defaults).")
        self._tau0(strict=True)
        if not t["regs"]: raise ValueError("Add at least one region.")
        if not t["cons"]: raise ValueError("Add at least one contact.")
        if min(t["Nx"],t["Ny"],t["Nz"])<2: raise ValueError("Each base node count must be at least 2.")
        if not min(t["Lx"],t["Ly"],t["Lz"])>0: raise ValueError("The domain sizes Lx, Ly, Lz must be positive.")
        if not 20.<=T<=2000.: raise ValueError(f"Temperature {T:g} K is outside 20 … 2000 K.")
        bad=[r.label for r in t["regs"] if not (r.x1>r.x0 and r.y1>r.y0 and r.z1>r.z0)]
        if bad: raise ValueError("Region(s) with an end coordinate not above the start: "+", ".join(bad))
        if len(t["cons"])>MAX_CONTACTS: raise ValueError(f"At most {MAX_CONTACTS} contacts.")
        models=self._models()
        graded=bool(self.use_graded.get())
        info=[f"'{c.label}': "+contact_core(c,t["regs"],(t["Lx"],t["Ly"],t["Lz"]),models)[2]
              for c in t["cons"] if getattr(c,"metal","")]
        if info and info!=getattr(self,"_metal_logged",None):
            self._metal_logged=info
            for x in info: self._log("Contact "+x)
        N=self._node_count()
        if N:
            if N*BYTES_PER_NODE>0.8*_ram_bytes():
                raise ValueError(f"The mesh has {N:,} nodes (~{N*BYTES_PER_NODE/1024**3:.1f} GB) - more than "
                                 f"this computer's free memory allows. Reduce Nx/Ny/Nz.")
        return t,solver,models,T,graded

    def _models(self):
        return ((M_BGN if self.m_bgn.get() else 0)|(M_TAU if self.m_tau.get() else 0)|
                (M_VSAT if self.m_vsat.get() else 0)|(M_QC if self.m_qc.get() else 0)|
                (M_LB if self.m_lb.get() else 0)|(M_TUN if self.m_tun.get() else 0)|
                (M_TE if self.m_te.get() else 0)|(M_SBT if self.m_sbt.get() else 0)|
                (M_II if self.m_ii.get() else 0))
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
                txt=getattr(job,"tag","")+f"point {min(job.k+1,job.N)}/{job.N}"+(f" (+{ex} refined)" if ex else "")
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
        amr=None
        if self.amr_on.get() and t.get("mesh")!="tri":
            try: amr=(max(1,int(self.amr_passes.get())),max(1.,float(self.amr_budget.get())),
                      max(0.05,float(self.amr_kappa.get())))
            except Exception:
                messagebox.showerror("Cannot run","Adaptive mesh refinement: passes, max nodes and threshold must be numbers."); return
        self.sv.set(f"Solving {self.tname} ({nn:,} nodes{', prisms' if t.get('mesh')=='tri' else ''}"
                    f"{', adaptive mesh' if amr else ''})…")
        def work(job):
            fresh=self._prepare(job,t,solver,models,T,graded,sig,reuse)
            info=[]
            if amr:
                base=dict(t); base["amr"]=None
                N0=int(np.prod([len(a) for a in make_grid(base,graded)]))
                ok,Vx,info=amr_run(t,cons,flt,job,vapp,solver,models,T,graded,passes=amr[0],
                                   max_nodes=int(amr[1]*N0),kappa=amr[2])
                if info:                                  # the C state is the refined device now
                    self._c_sig=device_signature(t,solver,models,T,graded); self._grid=tuple(make_grid(t,graded))
                    fresh=True
            else:
                ok,Vx=run_point(cons,flt,job,vapp)
            self._c_ok=True
            s=pull()
            if s is not None:
                s["dev"]=dict(regs=t["regs"],cons=cons,sheets=t.get("sheets",[]),L=(t["Lx"],t["Ly"],t["Lz"]),
                              defects=t.get("defects") or [])
                s["models"]=models
            return dict(ok=ok,Vx=Vx,s=s,fresh=fresh,dt=time.time()-job.t0,flt=flt,cons=cons,vapp=vapp,
                        its=api_stat(3),stopped=job.stop,notes=list(_DEFECT_NOTES) if t.get("defects") else [],
                        amr=copy.deepcopy(t.get("amr")) if amr else None,amr_info=info,amr_on=bool(amr))
        self._start(job,work,self._run_done)

    def _defect_notes(self,r):
        """Log what the defect mapping reported (defects that found no mesh node...)."""
        notes=r.get("notes") or []
        if r.get("fresh"): self._def_summary(notes)
        for n in notes if r.get("fresh") else []: self._log("Defect: "+n)
    def _run_done(self,r):
        s=r["s"]; self.sol=s; self._defect_notes(r)
        if r.get("amr_info"):
            self.amr_nodes=r.get("amr") or {}; self._amr_label(r["amr_info"])
            self._log("Adaptive mesh: "+"; ".join(f"pass {i+1} +{a} lines → {b:,} nodes{'' if c else ' (not converged)'}"
                                                  for i,(a,b,c) in enumerate(r["amr_info"])))
            self._mesh_update(redraw=False); self._ddev()
        elif r.get("amr_on") and r["ok"]:
            self._log("Adaptive mesh: the solution is resolved at this bias - no refinement needed")
        elif self.amr_on.get() and self.mesh_type.get()=="tri":
            self._log("Adaptive mesh refinement works on the rectangular mesh only - skipped")
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
                                            okey=okey,bv=bv,Vs=Vs,Ve=Ve,wl=self.tmeta.get("wl"),models=models,T=T)
        self.sv.set(f"Sweeping {tkey}: {Vs:g} → {Ve:g} V, {Np} points"+(f", {fkey} floating" if flt else "")+"…")
        self._log(f"Sweep {self.tname}: {tkey} {Vs:g}→{Ve:g} V ({Np} pts)"+(f", floating {fkey}" if flt else ""))
        refine=bool(self.iv_refine.get()); reuse=bool(self.reuse_sol.get())
        if self.iv_live.get(): self.nb.select(self.ti)
        def work(job):
            fresh=self._prepare(job,t,solver,models,T,graded,sig,reuse)
            R=run_sweep(cons,tgt,Vs,Ve,Np,flt,job,on_point=lambda k,rec:self._q.put(("point",rec)),
                        vapp=vapp,refine=refine,max_fail=2 if models&M_II else 0)
            self._c_ok=True
            A=np.array([api_conA(c) for c in range(len(cons))])
            s=pull()
            if s is not None:
                s["dev"]=dict(regs=t["regs"],cons=cons,sheets=t.get("sheets",[]),L=(t["Lx"],t["Ly"],t["Lz"]),
                              defects=t.get("defects") or [])
                s["models"]=models
            return dict(R=R,A=A,s=s,fresh=fresh,dt=time.time()-job.t0,stopped=job.stop,
                        cons=cons,vapp=vapp,tgt=tgt,flt=flt,notes=list(_DEFECT_NOTES) if t.get("defects") else [])
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
        R=r["R"]; m=self._live_meta; cons=r["cons"]; self._defect_notes(r)
        self.ivd=dict(R=R,A=r["A"],cons=cons,groups=con_groups(cons),tkey=m["tkey"],fkey=m["fkey"],
                      okey=m["okey"],bv=m["bv"],labels=[c.label for c in cons],bcs=[c.bc for c in cons],
                      vapp=r["vapp"],tname=self.tname,wl=m.get("wl"),models=m.get("models"),
                      Vs=m.get("Vs"),Ve=m.get("Ve"),T=m.get("T"))
        n=len(R["V"]); nc=int(np.sum(R["conv"])) if n else 0
        state="STOPPED" if r["stopped"] else "done"
        sf=R.get("stop_fail",float("nan"))
        if sf==sf:
            self._log(f"Sweep ended early at V = {sf:g} V: the last points did not converge"
                      +(" (avalanche breakdown)" if (m.get("models") or 0)&M_II else ""))
        self.cv2.set(f"sweep {nc}/{n} converged · {r['dt']:.1f} s")
        self.sv.set(f"Sweep {state}: {nc}/{n} points converged in {r['dt']:.1f} s (fields show the last point)")
        self._log(f"Sweep {state}: {nc}/{n} converged, {r['dt']:.1f} s")
        if r["s"] is not None:
            self.sol=r["s"]; self.ivd["s"]=r["s"]
            self._fill_results(r["s"],cons,r["vapp"],r["flt"],None,note="last sweep point")
            self._new_solution_views()
        self._plot_iv(); self.nb.select(self.ti)

    # ─── PINNED SWEEPS ─────────────────────────────────────────
    PIN_COLS=["#f472b6","#38bdf8","#facc15","#a78bfa","#fb923c","#4ade80","#f87171","#22d3ee"]
    PIN_FIELDS_MAX=600_000            # nodes: larger solutions are pinned without their fields
    def _pin_record(self,d,name,s=None,color=None):
        """A pinned copy of a sweep record (arrays and contacts copied)."""
        keep=dict(R={k:np.array(v,copy=True) for k,v in d["R"].items()},
                  A=None if d.get("A") is None else np.array(d["A"],copy=True),
                  cons=copy.deepcopy(d.get("cons")),groups=copy.deepcopy(d["groups"]),
                  tkey=d["tkey"],fkey=d.get("fkey"),okey=d.get("okey"),bv=d.get("bv",0.),wl=d.get("wl"),
                  models=d.get("models"),tname=d.get("tname",""),labels=list(d.get("labels") or []),
                  bcs=list(d.get("bcs") or [c.bc for c in d["cons"]]),vapp=list(d.get("vapp") or []),
                  Vs=d.get("Vs"),Ve=d.get("Ve"),T=d.get("T"))
        if s is not None and s.get("nodes",0)>self.PIN_FIELDS_MAX: s=None
        if color is None:
            used={p["color"] for p in self.pins}
            color=next((c for c in self.PIN_COLS if c not in used),self.PIN_COLS[len(self.pins)%len(self.PIN_COLS)])
        return dict(name=name,color=color,show=True,d=keep,s=s,when=time.strftime("%H:%M:%S"))
    def _unique_pin_name(self,name):
        have={p["name"] for p in self.pins}
        if name not in have: return name
        k=2
        while f"{name} ({k})" in have: k+=1
        return f"{name} ({k})"
    def _pin_sweep(self):
        d=self.ivd
        if d is None or not len(d["R"]["V"]):
            messagebox.showinfo("Pin sweep","Run a sweep first (↔ Sweep, F6): the last sweep is what gets pinned."); return
        name=self._unique_pin_name(d.get("tname") or self.tname or "sweep")
        self.pins.append(self._pin_record(d,name,s=d.get("s")))
        self._log(f"Pinned sweep '{name}'"); self.sv.set(f"Pinned '{name}' - the next sweeps are drawn over it")
        self._plot_iv(); self._pins_fill()
    def _pins_win(self):
        top=self._pin_top
        if top is not None and top.winfo_exists():
            top.deiconify(); top.lift(); self._pins_fill(); return
        top=tk.Toplevel(self); top.title("Pinned sweeps"); top.configure(bg=C["panel"]); top.transient(self)
        self._pin_top=top
        try: top.geometry(f"+{self.winfo_rootx()+120}+{self.winfo_rooty()+110}")
        except Exception: pass
        tk.Label(top,text="Pinned sweeps are drawn over the last sweep on the I-V tab. Click 'Shown' to hide or show one.",
                 bg=C["panel"],fg=C["sub"],font=self.F["small"],anchor="w").pack(fill=tk.X,padx=10,pady=(8,2))
        lf=tk.Frame(top,bg=C["panel"]); lf.pack(fill=tk.BOTH,expand=True,padx=8,pady=2)
        cols=(("show",56,"Shown"),("name",200,"Name"),("dev",190,"Device"),("sweep",210,"Sweep"),
              ("when",70,"Pinned"),("fields",56,"Fields"))
        tv=ttk.Treeview(lf,columns=[c for c,_,_ in cols],show="headings",height=10,selectmode="browse")
        for c,w,t in cols:
            tv.heading(c,text=t); tv.column(c,width=w,anchor="center" if c in ("show","when","fields") else "w",
                                            stretch=c in ("name","dev","sweep"))
        vs=ttk.Scrollbar(lf,orient="vertical",command=tv.yview); tv.configure(yscrollcommand=vs.set)
        tv.pack(side=tk.LEFT,fill=tk.BOTH,expand=True); vs.pack(side=tk.LEFT,fill=tk.Y)
        tv.bind("<ButtonRelease-1>",self._pins_click); tv.bind("<Double-1>",lambda e:self._pins_rename())
        self._pin_tv=tv
        bf=tk.Frame(top,bg=C["panel"]); bf.pack(fill=tk.X,padx=8,pady=(4,8))
        for txt,fn in (("Rename…",self._pins_rename),("Colour…",self._pins_colour),("Show fields",self._pins_fields),
                       ("Export CSV…",self._pins_csv),("Remove",self._pins_remove),("Remove all",self._pins_clear)):
            self._btn(bf,txt,fn,bg=C["border"]).pack(side=tk.LEFT,padx=3)
        self._pins_fill()
    def _pins_fill(self):
        self._ptab_fill()
        top=self._pin_top
        if top is None or not top.winfo_exists(): return
        tv=self._pin_tv; sel=tv.selection()
        for it in tv.get_children(): tv.delete(it)
        for k,p in enumerate(self.pins):
            d=p["d"]; R=d["R"]; V=R["V"]
            rng=f"{d['tkey']} {V.min():.3g}…{V.max():.3g} V ({int(np.sum(R['conv']))}/{len(V)} conv.)" if len(V) else d["tkey"]
            tag="p"+p["color"].lstrip("#"); tv.tag_configure(tag,foreground=p["color"])
            tv.insert("","end",iid=str(k),values=("✓" if p.get("show",True) else "–",p["name"],d.get("tname",""),rng,
                                                  p["when"],"yes" if p.get("s") is not None else "-"),tags=(tag,))
        if sel and tv.exists(sel[0]): tv.selection_set(sel[0])
        self.iv_pin_lbl.config(text=f"{len(self.pins)} pinned" if self.pins else "")
    def _pin_sel(self):
        tv=getattr(self,"_pin_tv",None)
        if tv is None or not tv.selection(): return None
        k=int(tv.selection()[0]); return k if k<len(self.pins) else None
    def _pins_click(self,e):
        tv=self._pin_tv; it=tv.identify_row(e.y); col=tv.identify_column(e.x)
        if not it or col!="#1": return
        p=self.pins[int(it)]; p["show"]=not p.get("show",True)
        self._pins_fill(); self._replot_iv()
    def _pins_rename(self):
        k=self._pin_sel()
        if k is None: return
        from tkinter import simpledialog
        v=simpledialog.askstring("Rename pinned sweep","Name:",initialvalue=self.pins[k]["name"],parent=self._pin_top)
        if v and v.strip(): self.pins[k]["name"]=v.strip(); self._pins_fill(); self._replot_iv()
    def _pins_colour(self):
        k=self._pin_sel()
        if k is None: return
        from tkinter import colorchooser
        c=colorchooser.askcolor(color=self.pins[k]["color"],parent=self._pin_top,title="Curve colour")
        if c and c[1]: self.pins[k]["color"]=c[1]; self._pins_fill(); self._replot_iv()
    def _pins_fields(self):
        k=self._pin_sel()
        if k is None: messagebox.showinfo("Pinned sweeps","Select a pinned sweep first.",parent=self._pin_top); return
        p=self.pins[k]
        if p.get("s") is None:
            messagebox.showinfo("Pinned sweeps","This pin was stored without its fields (an A/B run with 'keep the "
                                "fields' off, or a very large mesh).",parent=self._pin_top); return
        d=p["d"]; self.sol=p["s"]
        self._fill_results(p["s"],d["cons"],d.get("vapp") or None,None,None,note=f"pinned '{p['name']}', last sweep point")
        self._fresh_cuts=False; self._new_solution_views()
        self.sv.set(f"Fields of '{p['name']}' (its last sweep point) in the solution tabs")
    def _pins_remove(self):
        k=self._pin_sel()
        if k is None: return
        del self.pins[k]; self._pins_fill(); self._replot_iv()
        if not self.pins and self.ivd is None: self._plot_iv()
    def _pins_clear(self):
        if not self.pins: return
        if not messagebox.askyesno("Pinned sweeps",f"Remove all {len(self.pins)} pinned sweeps?",parent=self._pin_top): return
        self.pins=[]; self._pins_fill()
        if self.ivd is not None: self._plot_iv()
        else: self._placeholder(self.fiv,self.cviv,"No sweep yet - set it up in the Sweep page and press ↔ Sweep (F6)")
        self.iv_pin_lbl.config(text="")
    def _pins_csv(self):
        if not self.pins: return
        p=filedialog.asksaveasfilename(defaultextension=".csv",filetypes=[("CSV","*.csv")],parent=self._pin_top,
                                       initialfile="pinned_sweeps.csv")
        if not p: return
        n=max(len(q["d"]["R"]["V"]) for q in self.pins); cols=[]; hdr=[]
        for q in self.pins:
            P=self._iv_prep(q["d"]); V=P["V"]; I=P["gI"][P["okey"]]
            pad=lambda a:np.concatenate([np.asarray(a,float),np.full(n-len(a),np.nan)])
            nm=q["name"].replace(",",";")
            cols+= [pad(V),pad(I),pad(P["cv"].astype(float))]
            hdr+=[f"{nm}: V_{P['tkey']}_V",f"{nm}: I_{P['okey']}_A",f"{nm}: converged"]
        np.savetxt(p,np.column_stack(cols),delimiter=",",fmt="%.8e",comments="",header=",".join(hdr))
        self.sv.set(f"Pinned sweeps → {p}"); self._log(f"Pinned sweeps exported to {p}")

    # ─── A/B COMPARE ───────────────────────────────────────────
    MODEL_NAMES={M_BGN:"BGN",M_TAU:"τ(N)",M_VSAT:"v_sat",M_QC:"MLDA",M_LB:"Lombardi",M_TUN:"tunnelling",
                 M_TE:"TE",M_SBT:"barrier tunnelling",M_II:"avalanche"}
    def _ab_win(self):
        top=self._ab_top
        if top is not None and top.winfo_exists():
            top.deiconify(); top.lift(); self._ab_fill(); return
        top=tk.Toplevel(self); top.title("Compare designs A and B"); top.configure(bg=C["panel"]); top.transient(self)
        self._ab_top=top
        try: top.geometry(f"+{self.winfo_rootx()+100}+{self.winfo_rooty()+80}")
        except Exception: pass
        tk.Label(top,text="Capture two designs, then run the Sweep page's sweep on both - their curves come out side "
                          "by side on the I-V tab (each one pinned). A design is everything in the editor: structure, "
                          "mesh, contact voltages, physics models and temperature, as they are when you press Capture. "
                          "Load a template (or edit the device) and capture A; change it and capture B; Run.",
                 bg=C["panel"],fg=C["sub"],font=self.F["small"],justify="left",anchor="w",wraplength=600
                 ).pack(fill=tk.X,padx=10,pady=(8,4))
        row=tk.Frame(top,bg=C["panel"]); row.pack(fill=tk.X,padx=8,pady=2)
        self.ab_name={}; self.ab_lbl={}
        for k,slot in enumerate(("A","B")):
            fr=tk.LabelFrame(row,text=f" Design {slot} ",bg=C["panel"],fg=C["accent"],font=self.F["uib"],bd=1,relief="groove")
            fr.grid(row=0,column=k,sticky="nsew",padx=4); row.columnconfigure(k,weight=1)
            v=tk.StringVar(value=(self.ab[slot] or {}).get("name","")); self.ab_name[slot]=v
            f1=tk.Frame(fr,bg=C["panel"]); f1.pack(fill=tk.X,padx=6,pady=(4,2))
            tk.Label(f1,text="Name",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT)
            e=tk.Entry(f1,textvariable=v,width=24,bg="#07111e",fg=C["text"],insertbackground=C["text"],relief="flat",
                       font=self.F["ui"]); e.pack(side=tk.LEFT,fill=tk.X,expand=True,padx=4)
            e.bind("<FocusOut>",lambda _e,sl=slot:self._ab_rename(sl)); e.bind("<Return>",lambda _e,sl=slot:self._ab_rename(sl))
            l=tk.Label(fr,text="",bg=C["panel"],fg=C["text"],font=self.F["small"],justify="left",anchor="w",wraplength=280)
            l.pack(fill=tk.X,padx=6,pady=2); self.ab_lbl[slot]=l
            f2=tk.Frame(fr,bg=C["panel"]); f2.pack(fill=tk.X,padx=6,pady=(2,6))
            self._btn(f2,"Capture the editor's device",lambda sl=slot:self._ab_capture(sl)).pack(side=tk.LEFT,padx=2)
            self._btn(f2,"Load into the editor",lambda sl=slot:self._ab_load(sl),bg=C["border"]).pack(side=tk.LEFT,padx=2)
        vf=tk.Frame(top,bg=C["panel"]); vf.pack(fill=tk.X,padx=10,pady=4); self._ab_vf=vf
        tk.Label(vf,text="Variants of the loaded template",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT)
        self.ab_var=tk.StringVar()
        self.ab_vcb=ttk.Combobox(vf,textvariable=self.ab_var,values=[],state="readonly",width=34,font=self.F["ui"])
        self.ab_vcb.pack(side=tk.LEFT,padx=4)
        for slot in ("A","B"):
            self._btn(vf,f"→ {slot}",lambda sl=slot:self._ab_variant(sl),bg=C["border"]).pack(side=tk.LEFT,padx=2)
        self.ab_sw=tk.Label(top,text="",bg=C["panel"],fg=C["yellow"],font=self.F["small"],justify="left",anchor="w",
                            wraplength=600)
        self.ab_sw.pack(fill=tk.X,padx=10,pady=2)
        bf=tk.Frame(top,bg=C["panel"]); bf.pack(fill=tk.X,padx=8,pady=(4,8))
        self._check(bf,"Keep the final fields of both (Pins… → Show fields)",self.ab_keep,lambda:None)
        self._btn(bf,"Close",top.destroy,bg=C["border"]).pack(side=tk.RIGHT,padx=3)
        self.ab_run_btn=self._btn(bf,"▶ Run A and B",self._ab_run,bg=C["green"],fg="#03130a")
        self.ab_run_btn.pack(side=tk.RIGHT,padx=3)
        self._ab_fill()
    def _ab_fill(self):
        top=self._ab_top
        if top is None or not top.winfo_exists(): return
        mn=self.MODEL_NAMES
        for slot in ("A","B"):
            a=self.ab[slot]
            if a is None: self.ab_lbl[slot].config(text="(empty - press Capture)",fg=C["sub"]); continue
            t,solver,models,T,graded=a["snap"]
            ml=", ".join(v for b,v in mn.items() if models&b) or "none"
            vv=", ".join(f"{c.label} {c.V:g} V" for c in t["cons"][:6])+(" …" if len(t["cons"])>6 else "")
            fk=a.get("fkey","(none)")
            self.ab_lbl[slot].config(fg=C["text"],text=f"{a['tname']} · {a['nodes']:,} nodes · {T:g} K · captured {a['when']}\n"
                                                       f"models: {ml}\ncontacts: {vv}"+
                                                       (f"\nopen (floating, I = 0): {fk}" if fk!="(none)" else ""))
        base=self.tname.split(" [")[0]
        vars_=[lab for lab,_ in TEMPLATE_VARIANTS.get(base,[])]
        self.ab_vcb.configure(values=vars_,state="readonly" if vars_ else "disabled")
        if vars_ and self.ab_var.get() not in vars_: self.ab_var.set(vars_[0])
        if not vars_: self.ab_var.set("(this template has no variants)")
        tk_=self._gmap.get(self.iv_tgt.get()); fk=self._gmap.get(self.iv_flt.get(),"(none)"); ok=self._gmap.get(self.iv_out.get(),"auto")
        try: rng=f"{float(self.ivS.get()):g} → {float(self.ivE.get()):g} V, {int(self.ivN.get())} points"
        except Exception: rng="(check the Sweep page)"
        self.ab_sw.config(text=f"Sweep (from the Sweep page): {tk_ or '?'} {rng}, observing {ok}. The terminal names "
                               f"must exist in both designs; an open (floating) terminal is taken from each design as "
                               f"captured.")
    def _ab_rename(self,slot):
        a=self.ab.get(slot); v=self.ab_name[slot].get().strip()
        if a is not None and v: a["name"]=v
    def _ab_capture(self,slot,name=None):
        try: snap=self._snapshot()
        except ValueError as e: messagebox.showerror("Capture design",str(e),parent=self._ab_top); return False
        js=json.loads(json.dumps(self._device_json()))
        nm=(name or self.ab_name[slot].get().strip() or f"{slot}: {self.tname}")
        fk=self._gmap.get(self.iv_flt.get(),"(none)")          # an open terminal is part of the design
        self.ab[slot]=dict(snap=snap,js=js,tname=self.tname,name=nm,when=time.strftime("%H:%M:%S"),
                           nodes=self._node_count() or 0,fkey=fk or "(none)")
        self.ab_name[slot].set(nm); self._ab_fill()
        self._log(f"A/B: design {slot} = '{nm}' ({self.tname})")
        return True
    def _ab_load(self,slot):
        a=self.ab.get(slot)
        if a is None or self.busy: return
        js=copy.deepcopy(a["js"])
        try: check_device_json(js); parts=self._device_parts(js)
        except Exception as e: messagebox.showerror("Load design",str(e),parent=self._ab_top); return
        self._apply_device_json(js,*parts,src=a["tname"],label=f"A/B design {slot}")
        self.sv.set(f"Design {slot} ('{a['name']}') is in the editor"); self._ab_fill()
    def _ab_variant(self,slot):
        if self.busy: return
        base=self.tname.split(" [")[0]; lab=self.ab_var.get()
        fn=dict(TEMPLATE_VARIANTS.get(base,[])).get(lab)
        if fn is None: return
        try: t=fn()
        except Exception: messagebox.showerror("Template variant",traceback.format_exc(),parent=self._ab_top); return
        keep=(self.iv_tgt.get(),self.iv_flt.get(),self.iv_out.get(),self.ivS.get(),self.ivE.get(),self.ivN.get())
        self._apply_template(f"{base} [{lab}]",t,base=base)
        # the variants share their terminals: keep the sweep set up on the Sweep
        # page (an open, floating terminal belongs to the variant itself)
        if keep[0] in self._gmap:
            self.iv_tgt.set(keep[0])
            self.iv_out.set(keep[2] if keep[2]=="auto" or keep[2] in self._gmap else "auto")
            self.ivS.set(keep[3]); self.ivE.set(keep[4]); self.ivN.set(keep[5])
        self._ab_capture(slot,name=lab)
    def _ab_run(self):
        if self.busy: return
        A,B=self.ab["A"],self.ab["B"]
        if A is None or B is None:
            messagebox.showinfo("A/B compare","Capture both designs first.",parent=self._ab_top); return
        try: Vs=float(self.ivS.get()); Ve=float(self.ivE.get()); Np=int(self.ivN.get()); bv=float(self.bv_thr.get())
        except Exception:
            messagebox.showerror("A/B compare","Start, end, points and BV knee (Sweep page) must be numbers.",parent=self._ab_top); return
        if not (1<=Np<=2000): messagebox.showerror("A/B compare","Points must be 1…2000.",parent=self._ab_top); return
        tkey=self._gmap.get(self.iv_tgt.get())
        okey=self._gmap.get(self.iv_out.get(),"auto")
        if not tkey: messagebox.showerror("A/B compare","Choose the swept terminal in the Sweep page.",parent=self._ab_top); return
        plan=[]
        for slot,a in (("A",A),("B",B)):
            cons=a["snap"][0]["cons"]; tgt=group_of(cons,tkey)
            fkey=a.get("fkey","(none)")
            flt=group_of(cons,fkey) if fkey!="(none)" else []
            if not tgt or (fkey!="(none)" and not flt):
                messagebox.showerror("A/B compare",f"Design {slot} ('{a['name']}') has no terminal "
                                     f"'{tkey if not tgt else fkey}'.",parent=self._ab_top); return
            if set(flt)&set(tgt):
                messagebox.showerror("A/B compare",f"Design {slot}: the open terminal cannot be the swept one.",
                                     parent=self._ab_top); return
            plan.append((slot,a,tgt,flt))
        keepf=bool(self.ab_keep.get()); refine=bool(self.iv_refine.get())
        job=Job("sweep",Np); job.nthreads=self._nthreads()
        self.sv.set(f"A/B compare: sweeping {tkey} {Vs:g} → {Ve:g} V on both designs…")
        self._log(f"A/B compare: '{A['name']}' vs '{B['name']}', {tkey} {Vs:g}→{Ve:g} V ({Np} pts)")
        def work(job):
            out=[]
            for slot,a,tgt,flt in plan:
                if job.stop: break
                t,solver,models,T,graded=a["snap"]; cons=t["cons"]
                job.tag=f"design {slot} · "; job.extra=0
                sig=device_signature(t,solver,models,T,graded)
                self._prepare(job,t,solver,models,T,graded,sig,False)
                vapp=[c.V for c in cons]
                R=run_sweep(cons,tgt,Vs,Ve,Np,flt,job,vapp=vapp,refine=refine,max_fail=2 if models&M_II else 0)
                self._c_ok=True
                Aa=np.array([api_conA(c) for c in range(len(cons))])
                s=None
                if keepf:
                    s=pull()
                    if s is not None:
                        s["dev"]=dict(regs=t["regs"],cons=cons,sheets=t.get("sheets",[]),L=(t["Lx"],t["Ly"],t["Lz"]),
                                      defects=t.get("defects") or [])
                        s["models"]=models
                out.append(dict(slot=slot,name=a["name"],tname=a["tname"],R=R,A=Aa,s=s,cons=cons,vapp=vapp,
                                models=models,flt=flt,fkey=a.get("fkey") if flt else None,T=T))
            return dict(out=out,stopped=job.stop,dt=time.time()-job.t0,tkey=tkey,okey=okey,bv=bv)
        self._start(job,work,self._ab_done)
    def _ab_done(self,r):
        outs=r["out"]
        if not outs: self.sv.set("A/B compare stopped"); return
        for p in self.pins: p["show"]=False             # the comparison on its own (older pins: Pins…)
        cols={"A":"#38bdf8","B":"#f472b6"}; msg=[]
        for o in outs:
            cons=o["cons"]
            d=dict(R=o["R"],A=o["A"],cons=cons,groups=con_groups(cons),tkey=r["tkey"],fkey=o["fkey"],okey=r["okey"],
                   bv=r["bv"],labels=[c.label for c in cons],bcs=[c.bc for c in cons],vapp=o["vapp"],
                   tname=o["tname"],models=o["models"],T=o["T"])
            self.pins.append(self._pin_record(d,self._unique_pin_name(o["name"]),s=o["s"],color=cols[o["slot"]]))
            n=len(o["R"]["V"]); msg.append(f"{o['slot']} {int(np.sum(o['R']['conv'])) if n else 0}/{n}")
        last=outs[-1]
        if last["s"] is not None:
            self.sol=last["s"]
            self._fill_results(last["s"],last["cons"],last["vapp"],last["flt"],None,
                               note=f"A/B design {last['slot']} '{last['name']}', last sweep point")
            self._fresh_cuts=True; self._new_solution_views()
        self.ivd=None
        state="STOPPED" if r["stopped"] else "done"
        self.sv.set(f"A/B compare {state}: {', '.join(msg)} points converged in {r['dt']:.0f} s - both curves are pinned")
        self._log(f"A/B compare {state}: {', '.join(msg)} converged, {r['dt']:.1f} s")
        self._plot_iv(); self._pins_fill(); self.nb.select(self.ti)

    # ─── RESULTS TABLE ─────────────────────────────────────────
    def _fill_results(self,s,cons=None,vapp=None,flt=None,Vx=None,note=""):
        for it in self.res_tree.get_children(): self.res_tree.delete(it)
        if s is None:
            self.res_title.config(text="No solution yet - press ▶ Run (F5)"); return
        cons=cons or self.cons; flt=flt or []
        nconv=s.get("conv",True)
        raw=s.get("Iraw") or []
        kir=""
        if raw:
            big=max([abs(x) for x in raw]+[1e-300]); fl=sum(s.get("Ifl") or [0.])
            kir=(f" · Kirchhoff check of the contact fluxes {abs(sum(raw))/big:.1e}" if big>10.*fl else
                 " · currents below the resolution floor")
        self.res_title.config(text=f"{self.tname} - {'converged' if nconv else 'NOT converged'}"
                                   f"{' ('+note+')' if note else ''} · {_mesh_desc(s)}{kir}")
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
        vis=self._vis("td")
        if not vis: self._sec_axes=None; self._no_panels(fig,self.cvd); return
        W={"3d":1.0 if "sec" in vis else 2.2,"sec":1.55 if "3d" in vis else 2.6,"leg":0.62}
        ax=row_axes(fig,[(k,W[k],"3d" if k=="3d" else None) for k in vis],
                    wspace=0.08 if "3d" in vis else 0.05)
        a3=ax.get("3d"); aLg=ax.get("leg")
        # section hidden: draw it on an axes of a figure that is never shown (the
        # legend still lists what the plane cuts)
        a2=ax.get("sec") or Figure().add_subplot()
        if aLg is not None: aLg.axis("off")
        a2.set_facecolor(C["bg"])
        hide=self.lay_hide if self.lay_sec.get() else set()          # the Layers… choices
        ext=lambda r,ax:{"x":(r.x0,r.x1),"y":(r.y0,r.y1),"z":(r.z0,r.z1)}[ax]
        tol=1e-9*max(Lx,Ly,Lz); lab_pts=[]
        for r in self.regs:
            n0,n1=ext(r,N)
            if not (n0-tol<=pos<=n1+tol) or id(r) in hide: continue
            fc,al=self._reg_color(r)
            if self.lay_sec.get() and id(r) in self.lay_alpha: al=self.lay_alpha[id(r)]
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
            if not (n0-tol<=pos<=n1+tol) or id(c) in hide: hidden+=1; continue
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
            if getattr(c,"metal",""): bc=f"{c.metal}, {bc}"
            hand.append(Line2D([0],[0],color=col,lw=4,label=f"{k}: {c.label} ({bc}, {c.V:g} V)"+(f" [{c.net}]" if c.net else "")))
        for sh in self.sheets:
            oth={"x":("y","z"),"y":("x","z"),"z":("x","y")}[sh.axis]
            rg={sh.axis:(sh.pos,sh.pos),oth[0]:(sh.a0,sh.a1),oth[1]:(sh.b0,sh.b1)}
            if N!=sh.axis and rg[N][0]-tol<=pos<=rg[N][1]+tol:
                a2.plot(rg[H],rg[Vv],color="#f0abfc",lw=2,ls="--",zorder=5)
                hand.append(Line2D([0],[0],color="#f0abfc",lw=2,ls="--",label=f"σ {sh.sigma:+.1e} q/cm²"))
        if self.defects and (self.lay_def.get() or not self.lay_sec.get()):
            self._draw_defects_2d(a2,self.defects,H,Vv,N,pos,L3)
            hand.extend(self._defect_legend(self.defects))
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
        if hand and aLg is not None:
            if hidden:
                k=next((i for i,h in enumerate(hand) if h.get_label().startswith("Defects (")),len(hand))
                hand.insert(k,Line2D([0],[0],color="none",label=f"(+{hidden} contact(s) off this plane or hidden)"))
            aLg.legend(handles=hand,loc="upper left",fontsize=7,borderaxespad=0.,framealpha=0.9,
                       title="Contacts & defects" if self.defects else "Contacts",title_fontsize=8)
        self._sec_axes=(a2,H,Vv,N,pos) if "sec" in vis else None
        self._hover_txt=a2.text(0.01,0.01,"",transform=a2.transAxes,fontsize=7,color=C["yellow"],
                                va="bottom",ha="left",zorder=10,
                                bbox=dict(boxstyle="round,pad=0.25",facecolor=C["panel"],edgecolor=C["border"],alpha=0.85))
        self._hover_txt.set_visible(False)
        if not getattr(self,"_hover_cid",None):
            self._hover_cid=self.cvd.mpl_connect("motion_notify_event",self._struct_hover)
        if a3 is not None:
            a3.set_facecolor(C["bg"])
            self._draw_struct3(a3,self.regs,self.cons,L3)
            if self.lay_def.get(): self._draw_defects_3d(a3,self.defects,L3)
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

    # ─── LAYERS (show / hide parts of the Structure view) ───────
    def _lay_part(self,kind,k,obj=None):
        """(hidden, own opacity or None) of region / contact k as the Layers
        window set it.  obj: a solution's copy of the part - matched to the
        editor's part by position and label."""
        lst=self.regs if kind=="r" else self.cons
        if k>=len(lst) or (obj is not None and lst[k].label!=obj.label): return False,None
        o=lst[k]
        return id(o) in self.lay_hide,self.lay_alpha.get(id(o))
    def _reg_alpha(self,r):
        """A region's opacity in 3-D without its own setting (the global slider)."""
        ra=max(0.,min(1.,float(self.lay_ra.get()))); al=self._reg_color(r)[1]
        return min(al,ra) if ra<=0.35 else ra*min(1.,al/0.35)
    def _con_alpha(self,c):
        ca=max(0.,min(1.,float(self.lay_ca.get())))
        return 0.45*ca/0.85 if (c.face==6 and getattr(c,"shape",None) is not None) else ca
    def _draw_struct3(self,a3,regs,cons,L3,overlay=None,editor=True):
        """Regions and contacts in a 3-D axes (plot coordinates x, z, -y) as the
        Layers window shows them: hidden parts left out, each part's own
        opacity or the global sliders.  overlay: an opacity for every part
        without its own (the structure under the Current tab's arrows)."""
        Lx,Ly,Lz=L3
        ext=lambda r,ax:{"x":(r.x0,r.x1),"y":(r.y0,r.y1),"z":(r.z0,r.z1)}[ax]
        for k,r in enumerate(regs):
            hid,own=self._lay_part("r",k,None if editor else r)
            if hid or (r.mat in INSULATORS and not self.lay_ins3d.get()): continue
            fc,_=self._reg_color(r)
            a=own if own is not None else (overlay if overlay is not None else self._reg_alpha(r))
            a=max(0.,min(1.,a))
            faces=(_prism_faces(r.shape,{q:ext(r,q) for q in "xyz"}) if r.shape is not None
                   else self._box_faces(r.x0,r.x1,r.y0,r.y1,r.z0,r.z1))
            a3.add_collection3d(Poly3DCollection(faces,alpha=a,facecolor=fc,edgecolor=C["border"],
                                                 linewidth=0.25 if a>0.02 else 0.4))
        for k,c in enumerate(cons):
            hid,own=self._lay_part("c",k,None if editor else c)
            if hid: continue
            e=con_extent(c,Lx,Ly,Lz); col=self._CC[k%len(self._CC)]
            if c.face==6 and getattr(c,"shape",None) is not None: faces=_prism_faces(c.shape,e,caps=not c.shape.hole)
            else: faces=self._box_faces(e["x"][0],e["x"][1],e["y"][0],e["y"][1],e["z"][0],e["z"][1])
            a=own if own is not None else (min(1.,1.6*overlay) if overlay is not None else self._con_alpha(c))
            a=max(0.,min(1.,a))
            a3.add_collection3d(Poly3DCollection(faces,alpha=a,facecolor=col,edgecolor=col,linewidth=0.4 if a>0.02 else 0.8))
    def _layers_win(self):
        top=self._lay_top
        if top is not None and top.winfo_exists():
            top.deiconify(); top.lift(); self._lay_fill(); return
        top=tk.Toplevel(self); top.title("Layers - Structure and Current views"); top.configure(bg=C["panel"])
        top.transient(self); self._lay_top=top
        try: top.geometry(f"+{self.winfo_rootx()+80}+{self.winfo_rooty()+90}")
        except Exception: pass
        tk.Label(top,text="Click 'Shown' to hide or show a part. Select parts (Ctrl / Shift for several) to give them "
                          "their own opacity.",bg=C["panel"],fg=C["sub"],font=self.F["small"],anchor="w",justify="left",
                 wraplength=640).pack(fill=tk.X,padx=10,pady=(8,2))
        lf=tk.Frame(top,bg=C["panel"]); lf.pack(fill=tk.BOTH,expand=True,padx=8,pady=2)
        cols=(("show",54,"Shown"),("part",210,"Part"),("kind",80,"Kind"),("opac",70,"Opacity"),("info",230,"Details"))
        tv=ttk.Treeview(lf,columns=[c for c,_,_ in cols],show="headings",height=14,selectmode="extended")
        for c,w,t in cols:
            tv.heading(c,text=t); tv.column(c,width=w,anchor="center" if c in ("show","opac") else "w",
                                            stretch=c in ("part","info"))
        vs=ttk.Scrollbar(lf,orient="vertical",command=tv.yview); tv.configure(yscrollcommand=vs.set)
        tv.pack(side=tk.LEFT,fill=tk.BOTH,expand=True); vs.pack(side=tk.LEFT,fill=tk.Y)
        tv.bind("<ButtonRelease-1>",self._lay_click); tv.bind("<space>",self._lay_click)
        tv.bind("<<TreeviewSelect>>",self._lay_selected)
        self._lay_tv=tv
        sf=tk.Frame(top,bg=C["panel"]); sf.pack(fill=tk.X,padx=10,pady=(4,2)); sf.columnconfigure(1,weight=1)
        tk.Label(sf,text="Opacity of the selected parts",bg=C["panel"],fg=C["accent"],font=self.F["ui"]
                 ).grid(row=0,column=0,sticky="w",padx=(0,8))
        self._lay_sc=ttk.Scale(sf,from_=0.,to=1.,orient=tk.HORIZONTAL,variable=self.lay_sa,command=self._lay_sel_alpha)
        self._lay_sc.grid(row=0,column=1,sticky="ew")
        self._lay_sv=tk.Label(sf,text="",width=5,bg=C["panel"],fg=C["text"],font=self.F["ui"]); self._lay_sv.grid(row=0,column=2,padx=4)
        b=self._btn(sf,"Auto",self._lay_sel_auto,bg=C["border"]); b.grid(row=0,column=3,padx=2)
        Tip(b,"The selected parts follow the region and contact sliders again.",self.F["small"])
        bf=tk.Frame(top,bg=C["panel"]); bf.pack(fill=tk.X,padx=8,pady=4)
        for txt,fn in (("All visible",lambda:self._lay_set("all")),
                       ("Semiconductors only",lambda:self._lay_set("semi")),
                       ("Hide contacts",lambda:self._lay_set("nocon")),
                       ("See-through contacts",lambda:self._lay_set("glass"))):
            self._btn(bf,txt,fn,bg=C["border"]).pack(side=tk.LEFT,padx=3)
        of=tk.Frame(top,bg=C["panel"]); of.pack(fill=tk.X,padx=10,pady=(2,4)); of.columnconfigure(1,weight=1)
        redraw=lambda *_:self._debounce("_lay_after_id",self._lay_redraw,150)
        for r,(txt,var) in enumerate((("Region opacity (all)",self.lay_ra),("Contact opacity (all)",self.lay_ca))):
            tk.Label(of,text=txt,bg=C["panel"],fg=C["text"],font=self.F["ui"]).grid(row=r,column=0,sticky="w",padx=(0,8))
            ttk.Scale(of,from_=0.,to=1.,orient=tk.HORIZONTAL,variable=var,command=lambda *_:(redraw(),self._lay_fill_later())
                      ).grid(row=r,column=1,sticky="ew",pady=1)
        cf=tk.Frame(top,bg=C["panel"]); cf.pack(fill=tk.X,padx=6,pady=(0,8))
        self._check(cf,"Insulators in 3-D",self.lay_ins3d,self._lay_redraw)
        self._check(cf,"Defects",self.lay_def,self._lay_redraw)
        self._check(cf,"Apply to the cross-section too",self.lay_sec,self._lay_redraw)
        self._lay_fill()
    def _lay_redraw(self):
        """Redraw the views the Layers choices act on."""
        self._ddev()
        if self.sol is not None and (self.cur_st3.get() or self.cur_bg.get()!="|J|"): self._upd_curr_plot()
    def _lay_fill_later(self):
        self._debounce("_lay_fill_id",self._lay_fill,250)
    def _lay_fill(self):
        top=self._lay_top
        if top is None or not top.winfo_exists(): return
        tv=self._lay_tv; sel=tv.selection()
        for it in tv.get_children(): tv.delete(it)
        for k,r in enumerate(self.regs):
            fc,_=self._reg_color(r); tag="c"+fc.lstrip("#"); tv.tag_configure(tag,foreground=fc)
            ins=r.mat in INSULATORS
            info=r.mat+("" if ins else f" {r.dtype} {r.doping:.1e} cm⁻³")
            hid=id(r) in self.lay_hide or (ins and not self.lay_ins3d.get())
            a=self.lay_alpha.get(id(r)); op=f"{a:.2f}" if a is not None else f"auto {self._reg_alpha(r):.2f}"
            tv.insert("","end",iid=f"r{k}",values=("–" if hid else "✓",r.label,"insulator" if ins else "region",op,info),
                      tags=(tag,))
        for k,c in enumerate(self.cons):
            col=self._CC[k%len(self._CC)]; tag="c"+col.lstrip("#"); tv.tag_configure(tag,foreground=col)
            bc={0:"ohmic",1:"Schottky",2:"gate"}[c.bc]
            a=self.lay_alpha.get(id(c)); op=f"{a:.2f}" if a is not None else f"auto {self._con_alpha(c):.2f}"
            tv.insert("","end",iid=f"c{k}",values=("–" if id(c) in self.lay_hide else "✓",c.label,"contact",op,
                                                   f"{Con.FNAME[c.face]}, {bc}, {c.V:g} V"+(f" [{c.net}]" if c.net else "")),
                      tags=(tag,))
        if self.defects:
            tv.insert("","end",iid="d",values=("✓" if self.lay_def.get() else "–",f"{len(self.defects)} defect(s)","defects",
                                               "","the Defects page"))
        keep=[i for i in sel if tv.exists(i)]
        if keep: tv.selection_set(keep)
    def _lay_obj(self,it):
        if not it or it=="d": return None
        k=int(it[1:]); lst=self.regs if it[0]=="r" else self.cons
        return lst[k] if k<len(lst) else None
    def _lay_click(self,e=None):
        tv=self._lay_tv
        if e is not None and e.type!=tk.EventType.KeyPress:
            if tv.identify_column(e.x)!="#1": return            # elsewhere: just select
            items=[tv.identify_row(e.y)]
        else: items=list(tv.selection())
        changed=False
        for it in items:
            if not it: continue
            if it=="d": self.lay_def.set(not self.lay_def.get()); changed=True; continue
            o=self._lay_obj(it)
            if o is None: continue
            if it[0]=="r" and o.mat in INSULATORS and not self.lay_ins3d.get():
                self.lay_ins3d.set(True); self.lay_hide.discard(id(o))
            elif id(o) in self.lay_hide: self.lay_hide.discard(id(o))
            else: self.lay_hide.add(id(o))
            changed=True
        if changed: self._lay_fill(); self._lay_redraw()
    def _lay_selected(self,_=None):
        """Show the opacity of the first selected part on the slider."""
        objs=[self._lay_obj(i) for i in self._lay_tv.selection()]; objs=[o for o in objs if o is not None]
        if not objs: self._lay_sv.config(text=""); return
        o=objs[0]; a=self.lay_alpha.get(id(o))
        if a is None: a=self._reg_alpha(o) if isinstance(o,Reg) else self._con_alpha(o)
        self.lay_sa.set(round(a,2)); self._lay_sv.config(text=f"{a:.2f}")
    def _lay_sel_alpha(self,*_):
        objs=[self._lay_obj(i) for i in self._lay_tv.selection()]; objs=[o for o in objs if o is not None]
        if not objs: return
        a=max(0.,min(1.,float(self.lay_sa.get())))
        for o in objs: self.lay_alpha[id(o)]=a
        self._lay_sv.config(text=f"{a:.2f}")
        self._lay_fill_later(); self._debounce("_lay_after_id",self._lay_redraw,150)
    def _lay_sel_auto(self):
        for i in self._lay_tv.selection():
            o=self._lay_obj(i)
            if o is not None: self.lay_alpha.pop(id(o),None)
        self._lay_fill(); self._lay_selected(); self._lay_redraw()
    def _lay_set(self,what):
        if what=="all":
            self.lay_hide.clear(); self.lay_alpha.clear(); self.lay_def.set(True)
            self.lay_ca.set(0.85); self.lay_ra.set(0.35)
        elif what=="semi":
            self.lay_hide={id(o) for o in self.cons}|{id(r) for r in self.regs if r.mat in INSULATORS}
            self.lay_ins3d.set(False); self.lay_def.set(False)
        elif what=="nocon":
            self.lay_hide|={id(o) for o in self.cons}
        elif what=="glass":
            self.lay_hide-={id(o) for o in self.cons}; self.lay_ca.set(0.15)
            for o in self.cons: self.lay_alpha.pop(id(o),None)
        self._lay_fill(); self._lay_redraw()
    def _lay_prune(self):
        """Forget hidden parts and opacities of parts that no longer exist."""
        live={id(o) for o in self.regs}|{id(o) for o in self.cons}
        self.lay_hide&=live
        self.lay_alpha={k:v for k,v in self.lay_alpha.items() if k in live}
        self._lay_fill()

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
        hd=self._defect_at(P,H,Vv,N,pos,L) if self.defects else None
        if hd is not None: s+=f"   ◆ {hd.label}: {hd.strength_text()}"
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
        tr=ax.get_xaxis_transform(); k=0; last=[-1e30,-1e30]; row=0; used=0
        lo,hi=min(x0,x1),max(x0,x1)
        while k<n:
            j=k
            while j+1<n and own[j+1]==own[k]: j+=1
            idx=int(own[k])
            if 0<=idx<len(regs):
                r=regs[idx]; fc,_=self._reg_color(r); a,b=mids[k],mids[j+1]
                ax.axvspan(a,b,color=fc,alpha=0.12,lw=0,zorder=0)
                if k>0: ax.axvline(a,color=fc,lw=0.7,alpha=0.55,zorder=0.5)
                lab=r.label; tw,_=self._txt_px(ax,lab,6.)
                xc=min(max(0.5*(a+b),lo),hi)
                cpx=(xc-x0)*ppd; room=abs(b-a)*ppd
                if self.show_labels.get() and room>=0.5*tw:   # fits (may overhang its span a little)
                    # a strip ABOVE the plot area (under the title): never on top of curves or legends
                    for rr in (row,1-row):             # first free row (no overlap with its last label)
                        if cpx-0.5*tw>last[rr]+4:
                            # rows a fixed number of points above the axes: the title pad below matches
                            ax.annotate(lab,xy=(xc,1.0),xycoords=tr,xytext=(0,3+10*rr),textcoords="offset points",
                                        ha="center",va="bottom",fontsize=6,color=fc,annotation_clip=False,zorder=6,
                                        bbox=dict(boxstyle="round,pad=0.12",facecolor="#05090f",alpha=0.8,edgecolor="none"))
                            last[rr]=cpx+0.5*tw; row=1-rr; used=max(used,rr+1); break
            k=j+1
        if used:                                       # lift the title above the label rows
            ax.set_title(ax.get_title(),fontsize=9,fontweight="bold",pad=6+10*used)
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
    def _ov_plane(self,ax,hc,vc,own2,regs,cons=None,axes=None,pos=None,L=None,cuts=()):
        """Region boundaries (the cell faces between nodes of different
        regions: exactly what the solver sees), region labels and the
        contacts cut by the plane on a 2-D map.  own2: shape (len(vc), len(hc)).
        cuts: [('h', y) | ('v', x)] line cuts drawn on the map - labels keep
        clear of them.  Call after the axis limits (and any inversion) are set."""
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
            ax.add_collection(LineCollection(np.concatenate(segs),colors="white",linewidths=0.6,alpha=0.75,zorder=4))
        cpaths=[]                                   # contact areas: no labels on them
        if cons and axes and L:
            Hn,Vn,Nn=axes; tol=1e-9*max(L)
            for k,c in enumerate(cons):
                e=con_extent(c,*L); col=self._CC[k%len(self._CC)]
                if not (e[Nn][0]-tol<=pos<=e[Nn][1]+tol): continue
                h0,h1=e[Hn]; v0,v1=e[Vn]
                if c.face==6:
                    if getattr(c,"shape",None) is not None:
                        pgs=_section_polys(c.shape,e,Hn,Vn,Nn,pos) or []
                    else:
                        pgs=[(np.array([(h0,v0),(h1,v0),(h1,v1),(h0,v1)]),[])]
                    for pg in pgs:
                        # light hatch + crisp outline: the metal stays visible without hiding the map
                        ax.add_patch(_poly_patch(pg,facecolor="none",edgecolor=col,lw=0.,hatch="//",alpha=0.35,zorder=5))
                        p=_poly_patch(pg,facecolor="none",edgecolor=col,lw=1.0,alpha=0.9,zorder=6)
                        ax.add_patch(p); cpaths.append(p.get_path())
                elif abs(h1-h0)<tol or abs(v1-v0)<tol:
                    ax.plot([h0,h1],[v0,v1],color=col,lw=4,solid_capstyle="butt",zorder=6,clip_on=False)
        defs=((self.sol or {}).get("dev") or {}).get("defects")
        if defs and axes and L:
            self._draw_defects_2d(ax,defs,axes[0],axes[1],axes[2],pos,tuple(L),overlay=True)
        if self.show_labels.get():
            ax._ov_label_job=(hm,vm,own2,regs,cuts,cpaths)   # placed by _flush_labels once the layout is final

    def _flush_labels(self,fig):
        """Place the pending region labels of every map in fig (call after
        tight_layout/colorbars, just before drawing: placement is in pixels).
        The canvas re-places them whenever a map's size changes (window
        resize, first display of a tab) - see _labels_on_draw."""
        for ax in fig.axes:
            if getattr(ax,"_ov_label_job",None) is not None: self._place_labels(ax)

    def _place_labels(self,ax):
        for a in getattr(ax,"_ov_label_artists",None) or []:
            try: a.remove()
            except Exception: pass
        ax._ov_label_artists=[]
        try: ax.apply_aspect()
        except Exception: pass
        bb=ax.get_window_extent(); ax._ov_label_size=(bb.width,bb.height)
        try: ax._ov_label_artists=self._ov_labels(ax,*ax._ov_label_job) or []
        except Exception: traceback.print_exc()

    def _clear_3d_on_draw(self,ev):
        """draw_event: a 3-D view whose projected box, tick labels or axis
        labels run into a neighbouring plot (a long, flat device seen in
        perspective pokes out of its grid cell) is moved clear of it - into
        free room left/right/below/above - and shrunk only if there is not
        enough room (its labels keep their distance from the box, so shrinking
        alone does not help); then one more draw.  Its position is restored
        first whenever the canvas changes size."""
        fig=ev.canvas.figure; r=ev.renderer; moved=False
        FW,FH=fig.bbox.width,fig.bbox.height
        key=(round(FW),round(FH))
        a3s=[a for a in fig.axes if getattr(a,"name","")=="3d" and a.get_visible()]
        for a3 in a3s:
            if getattr(a3,"_clr_key",None)!=key:
                p0=getattr(a3,"_clr_pos0",None)
                if p0 is None: a3._clr_pos0=a3.get_position().frozen()
                elif a3._clr_key is not None: a3.set_position(p0); moved=True
                a3._clr_key=key; a3._clr_n=0
                if moved: continue
            if a3._clr_n>=8: continue
            try:
                tb=a3.get_tightbbox(r)
                obs=[o.get_tightbbox(r) for o in fig.axes if o is not a3 and o.get_visible()]
            except Exception: continue
            if tb is None: continue
            obs=[o for o in obs if o is not None and o.width>1 and o.height>1]
            hit=[o for o in obs if tb.overlaps(o)]
            if not hit: continue
            o=hit[0]; g=6.
            band_y=[b for b in obs if b.y1>tb.y0 and b.y0<tb.y1 and b is not o]
            band_x=[b for b in obs if b.x1>tb.x0 and b.x0<tb.x1 and b is not o]
            roomL=tb.x0-max([0.]+[b.x1 for b in band_y if b.x1<=tb.x0+0.5])
            roomR=min([FW]+[b.x0 for b in band_y if b.x0>=tb.x1-0.5])-tb.x1
            roomD=tb.y0-max([0.]+[b.y1 for b in band_x if b.y1<=tb.y0+0.5])
            roomU=min([FH]+[b.y0 for b in band_x if b.y0>=tb.y1-0.5])-tb.y1
            cands=[(tb.x1-o.x0+g,roomL,-1,0),(o.x1-tb.x0+g,roomR,1,0),
                   (tb.y1-o.y0+g,roomD,0,-1),(o.y1-tb.y0+g,roomU,0,1)]
            fit=[c for c in cands if 0<c[0]<=c[1]]
            if fit: need,room,sx,sy=min(fit,key=lambda c:c[0]); shift,f=need,1.
            else:   need,room,sx,sy=max(cands,key=lambda c:c[1]-c[0]); shift,f=max(room-1.,0.),0.85
            p=a3.get_position(); w,h=p.width*f,p.height*f
            x0=p.x0+sx*shift/FW+(p.width-w if sx>0 else 0.)+(0.5*(p.width-w) if sx==0 else 0.)
            y0=p.y0+sy*shift/FH+(p.height-h if sy>0 else 0.)+(0.5*(p.height-h) if sy==0 else 0.)
            a3.set_position([x0,y0,w,h]); a3._clr_n+=1; moved=True
        if moved: self._redraw_soon(ev.canvas)

    def _redraw_soon(self,cv):
        """One more draw once this one is over (TkAgg ignores draw_idle()
        from inside its own draw).  Requests are merged, and a canvas gets at
        most 6 chained re-draws in 3 s, so the draw handlers can never keep
        it redrawing (and the window busy) for ever."""
        t=time.monotonic(); n,t0=getattr(cv,"_chain",(0,0.))
        if t-t0>3.: n,t0=0,t
        if n>=6: return
        cv._chain=(n+1,t0)
        if getattr(cv,"_redraw_job",None) is None:
            def go():
                cv._redraw_job=None
                try: cv.draw_idle()
                except Exception: pass
            cv._redraw_job=self.after(20,go)

    def _labels_on_draw(self,ev):
        """draw_event: a map drawn at another size than its labels were placed
        for gets them placed again (then one more draw)."""
        redo=False
        for ax in ev.canvas.figure.axes:
            if getattr(ax,"_ov_label_job",None) is None: continue
            bb=ax.get_window_extent(); w0,h0=getattr(ax,"_ov_label_size",(0,0))
            if abs(bb.width-w0)>2 or abs(bb.height-h0)>2:
                self._place_labels(ax); redo=True
        if redo: self._redraw_soon(ev.canvas)

    def _ov_labels(self,ax,hm,vm,own2,regs,cuts,cpaths):
        """Region labels that stay readable: each label goes where its whole
        box fits inside its own region, clear of region boundaries, line cuts,
        contacts and the labels already placed (a pixel raster of the axes with
        a summed-area fit test), as central as possible; a label that fits
        nowhere is left out (hover over the map to read any region)."""
        bb=ax.get_window_extent(); W=int(bb.width); H=int(bb.height)
        if W<60 or H<40: return
        st=max(1,int(np.ceil(max(W,H)/480.)))       # raster step (screen px)
        Wr,Hr=max(W//st,8),max(H//st,8)
        x0,x1=ax.get_xlim(); y0,y1=ax.get_ylim()
        px=x0+(np.arange(Wr)+0.5)/Wr*(x1-x0); py=y0+(np.arange(Hr)+0.5)/Hr*(y1-y0)   # row 0 = bottom edge
        lo_h,hi_h=min(hm[0],hm[-1]),max(hm[0],hm[-1]); lo_v,hi_v=min(vm[0],vm[-1]),max(vm[0],vm[-1])
        ih=np.clip(np.searchsorted(hm,px,side="right")-1,0,len(hm)-2)
        jv=np.clip(np.searchsorted(vm,py,side="right")-1,0,len(vm)-2)
        R=np.asarray(own2)[jv[:,None],ih[None,:]]
        R=np.where(((py>=lo_v)&(py<=hi_v))[:,None]&((px>=lo_h)&(px<=hi_h))[None,:],R,-2)
        blk=np.zeros((Hr,Wr),bool)
        d=R[:,1:]!=R[:,:-1]; blk[:,1:]|=d; blk[:,:-1]|=d
        d=R[1:,:]!=R[:-1,:]; blk[1:,:]|=d; blk[:-1,:]|=d
        for kind,val in cuts or ():
            if val is None: continue
            if kind=="v":
                c=int((val-x0)/((x1-x0) or 1.)*Wr)
                if 0<=c<Wr: blk[:,max(0,c-2):c+3]=True
            else:
                r=int((val-y0)/((y1-y0) or 1.)*Hr)
                if 0<=r<Hr: blk[max(0,r-2):r+3,:]=True
        if cpaths:
            G=np.stack(np.meshgrid(px,py),-1).reshape(-1,2)
            for pth in cpaths:
                try: blk|=pth.contains_points(G).reshape(Hr,Wr)
                except Exception: pass
        ids,cnt=np.unique(R[R>=0],return_counts=True)
        arts=[]
        for idx in ids[np.argsort(-cnt)]:
            if idx>=len(regs): continue
            lab=regs[idx].label; own_free=(R==idx)&~blk
            if own_free.sum()<20: continue
            for fs,mg in ((6.5,3),(6.5,1),(5.5,1)):
                tw,th=self._txt_px(ax,lab,fs)
                bw=int(np.ceil((tw+4)/st))+2*mg; bh=int(np.ceil((th+3)/st))+2*mg
                if bw>=Wr or bh>=Hr: continue
                S=np.zeros((Hr+1,Wr+1),np.int32); S[1:,1:]=np.cumsum(np.cumsum(~own_free,0),1)
                tot=S[bh:,bw:]-S[:-bh,bw:]-S[bh:,:-bw]+S[:-bh,:-bw]
                rr,cc=np.nonzero(tot==0)
                if not len(rr): continue
                k=int(np.argmin((rr-rr.mean())**2+(cc-cc.mean())**2))
                r0,c0=rr[k],cc[k]
                xc=x0+(c0+bw/2)/Wr*(x1-x0); yc=y0+(r0+bh/2)/Hr*(y1-y0)
                arts.append(ax.text(xc,yc,lab,fontsize=fs,color="white",ha="center",va="center",zorder=7,clip_on=True,
                            bbox=dict(boxstyle="round,pad=0.18",facecolor="#05090f",alpha=0.72,edgecolor="none")))
                blk[max(0,r0-1):r0+bh+1,max(0,c0-1):c0+bw+1]=True
                break
        return arts

    def _reg_hover(self,ax,hc,vc,Z,lbl,axes,pos):
        """Remember a 2-D map for the hover readout (region + value under the cursor)."""
        if not hasattr(self,"_hmaps"): self._hmaps={}
        self._hmaps[ax]=dict(hc=np.asarray(hc,float),vc=np.asarray(vc,float),Z=np.asarray(Z),lbl=lbl,
                             axes=axes,pos=pos)

    def _map_hover(self,ev,label):
        """Status-line readout of the map under the cursor: coordinates, the
        plotted value, the region (material, doping) and any contact there."""
        m=getattr(self,"_hmaps",{}).get(ev.inaxes) if ev.inaxes is not None else None
        if m is None or ev.xdata is None or self.sol is None:
            if label.cget("text"): label.config(text="")
            return
        h,v=float(ev.xdata),float(ev.ydata); hc,vc,Z=m["hc"],m["vc"],m["Z"]
        i=int(np.argmin(np.abs(hc-h))); j=int(np.argmin(np.abs(vc-v)))
        try: val=float(Z[j,i])
        except Exception: val=float("nan")
        Hn,Vn,Nn=m["axes"]
        txt=f"{Hn} {h:.2f}  {Vn} {v:.2f} nm   {m['lbl']} = "+(f"{val:.4g}" if np.isfinite(val) else "—")
        dev=self.sol.get("dev") or {}; regs=dev.get("regs") or []; L=dev.get("L")
        if regs and L:
            P={Hn:np.array([h]),Vn:np.array([v]),Nn:np.array([m["pos"]])}
            try:
                ri=int(np.asarray(SM.region_index(regs,P["x"],P["y"],P["z"],L)).ravel()[0])
                if 0<=ri<len(regs):
                    r=regs[ri]
                    txt+=f"   │ {r.label}: {r.mat}"+("" if r.mat in INSULATORS else f" {r.dtype} {r.doping:.2e} cm⁻³")
                for k,c in enumerate(dev.get("cons") or []):
                    if c.face==6 and bool(np.asarray(SM.box_contact_mask(c,P["x"],P["y"],P["z"],L)).ravel()[0]):
                        txt+=f"   │ contact {c.label}"
            except Exception: pass
        label.config(text=txt)

    # ─── BANDS ─────────────────────────────────────────────────
    def _plot_band(self):
        s=self.sol; fig=self.fb; fig.clear()
        if s is None: self.cvb.draw_idle(); return
        vis=self._vis("tb")
        if not vis: self._no_panels(fig,self.cvb); return
        axs=row_axes(fig,[(k,1.6 if k=="band" else 1.,"3d" if k=="ribbon" else None) for k in vis])
        ax=axs.get("band") or Figure().add_subplot()       # hidden: drawn off screen
        ax3=axs.get("ribbon")
        ax.set_facecolor(C["bg"])
        if ax3 is not None: ax3.set_facecolor(C["bg"])
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
        if ax3 is not None:
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
        self._hmaps={k:v for k,v in getattr(self,"_hmaps",{}).items() if k.figure is not fig}
        self.cur_hover.config(text="")
        if s is None: self.cvc.draw_idle(); return
        vis=self._vis("tc")
        pkey=self._gmap.get(self.cur_prof.get()) if self.cur_prof.get()!="(none)" else None
        if pkey is None and "prof" in vis: vis.remove("prof")
        if not vis: self._no_panels(fig,self.cvc); return
        axs=grid_axes(fig,[(k,"3d" if k=="3d" else None) for k in vis])
        aJ,aR,aQ,aQ3,aP=(axs.get(k) for k in ("line","rec","map","3d","prof"))
        for a in axs.values(): a.set_facecolor(C["bg"])
        A,p1,p2=self.cur_axis.get(),self.cur_p1.get(),self.cur_p2.get()
        cx,Jx,(an1,av1),(an2,av2)=self._line_cut_at(s["Jx"],A,p1,p2)
        Jy=self._line_cut_at(s["Jy"],A,p1,p2)[1]; Jz=self._line_cut_at(s["Jz"],A,p1,p2)[1]
        Rt=self._line_cut_at(s["R"],A,p1,p2)[1]
        self.cur_lbl1.configure(text=f"{an1} = {av1:.1f} nm"); self.cur_lbl2.configure(text=f"{an2} = {av2:.1f} nm")
        xl=f"{'depth y' if A=='Y' else A.lower()} (nm)"
        own=self._owner(s); dev=s.get("dev") or {}; regs=dev.get("regs",[])
        ownl=self._line_cut_at(own,A,p1,p2)[1] if own is not None else None
        if aJ is not None:
            aJ.plot(cx,Jx,color="#34d399",lw=1.6,label="Jx"); aJ.plot(cx,Jy,color="#fbbf24",lw=1.6,label="Jy")
            aJ.plot(cx,Jz,color="#f87171",lw=1.6,label="Jz")
            aJ.plot(cx,np.sqrt(Jx**2+Jy**2+Jz**2),color=C["green"],lw=2.2,ls="--",label="|J|")
            aJ.axhline(0,color=C["border"],lw=0.5); aJ.legend(fontsize=7,framealpha=0.85)
            sax(aJ,xl,"J (A/cm²)",f"Current density along {A} ({an1} = {av1:.1f}, {an2} = {av2:.1f} nm)")
            self._ov_line(aJ,cx,ownl,regs)
        if aR is not None:
            if np.any(np.abs(np.nan_to_num(Rt))>1e-30):
                aR.plot(cx,Rt,color=C["red"],lw=2.,label="R net"); aR.legend(fontsize=7)
            else:
                aR.text(0.5,0.5,"R ≈ 0 (near equilibrium)",transform=aR.transAxes,ha="center",va="center",fontsize=10,color=C["sub"])
            sax(aR,xl,"R (cm⁻³s⁻¹)",f"Net recombination along {A}")
            self._ov_line(aR,cx,ownl,regs)
        if aP is not None: self._plot_prof(aP,s,pkey)
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
        n2,n3=self._arrow_n()
        bg=self.cur_bg.get(); sa=max(0.,min(1.,float(self.cur_sa.get())))
        if aQ is None: pass
        elif np.nanmax(M)>1e-30:
            X1,Y1=np.meshgrid(hx,vy)
            if bg!="structure":
                pc=aQ.pcolormesh(X1,Y1,np.ma.masked_invalid(M),cmap="viridis",shading="nearest")
                fig.colorbar(pc,ax=aQ,label="|J| (A/cm²)",shrink=0.85,pad=0.02)
            if bg!="|J|": self._struct_plane(aQ,hx,vy,own2,regs,sa)
            if self.cur_show_arrows2d.get():
                # n2 arrows across, rows spaced alike on screen; J interpolated to the
                # arrow points (so the count is not limited by the mesh)
                try:
                    W_,H_=fig.get_size_inches(); bb=aQ.get_position()
                    asp=(bb.height*H_)/max(bb.width*W_*0.85,1e-6)       # (the colour bar takes ~15 %)
                except Exception: asp=14/22.
                ny2=max(3,int(round(n2*min(max(asp,0.15),6.))))
                ex,ey=(hx[-1]-hx[0])/(2.*n2),(vy[-1]-vy[0])/(2.*ny2)   # cell centres
                gx=np.linspace(hx[0]+ex,hx[-1]-ex,n2); gy=np.linspace(vy[0]+ey,vy[-1]-ey,ny2)
                Us=_lin_axis(_lin_axis(U,hx,gx,1),vy,gy,0); Vs=_lin_axis(_lin_axis(Vq,hx,gx,1),vy,gy,0)
                ms=np.sqrt(Us**2+Vs**2)
                okm=np.isfinite(ms)&(ms>1e-30*np.nanmax(M)+1e-300)
                ln=np.log10(np.where(okm,ms,1e-300)); lo,hi=(np.nanmin(ln[okm]),np.nanmax(ln[okm])) if okm.any() else (0.,0.)
                sc=0.3+0.7*(ln-lo)/(hi-lo+1e-30); sc=np.where(okm,sc,0.)
                Gx,Gy=np.meshgrid(gx,gy)
                f=22./n2                     # denser arrows: shorter and thinner
                aQ.quiver(Gx,Gy,np.where(okm,Us/np.where(okm,ms,1),0)*sc,
                          np.where(okm,Vs/np.where(okm,ms,1),0)*sc*(-1 if invert else 1),
                          color="white",scale=25./f,width=min(0.006,max(0.0012,0.004*f**0.7)),
                          headwidth=4,headlength=5,alpha=0.9,pivot="mid")
        else:
            if bg!="|J|": self._struct_plane(aQ,hx,vy,own2,regs,sa)
            aQ.text(0.5,0.5,"J ≈ 0 (no current flow)",transform=aQ.transAxes,ha="center",va="center",fontsize=10,color=C["sub"])
        if aQ is not None:
            if cut[0]=="h": aQ.axhline(cut[1],color=C["yellow"],lw=0.9,ls="--",alpha=0.8)
            else: aQ.axvline(cut[1],color=C["yellow"],lw=0.9,ls="--",alpha=0.8)
            aQ.set_xlim(hx[0],hx[-1]); aQ.set_ylim(vy[0],vy[-1])
            if invert: aQ.invert_yaxis()
            self._ov_plane(aQ,hx,vy,own2,regs,dev.get("cons"),ovax,ovpos,dev.get("L"),[cut])
            sax(aQ,xlab,ylab,ttl+"  (dashed: the line cut)")
            self._reg_hover(aQ,hx,vy,M,"|J| (A/cm²)",ovax,ovpos)
        if aQ3 is not None and self.cur_st3.get() and regs:      # the structure, see-through, under the arrows
            L3d=dev.get("L") or (float(s["x"][-1]),float(s["y"][-1]),float(s["z"][-1]))
            self._draw_struct3(aQ3,regs,dev.get("cons") or [],L3d,overlay=sa,editor=False)
        if aQ3 is None: pass
        elif self.cur_show_arrows3d.get():
            x,y,z=s["x"],s["y"],s["z"]
            # a regular lattice of n3 x 0.75 n3 x 0.625 n3 arrows (x, depth y, z;
            # 8 x 6 x 5 by default), never more planes than mesh lines along an
            # axis (a 2-D cross-section keeps its few z planes), at most
            # ARROWS_3D_MAX arrows; J interpolated between the nodes
            Ls=np.array([x[-1]-x[0],y[-1]-y[0],z[-1]-z[0]],float); Lm=max(Ls.max(),1e-9)
            nn=(len(x),len(y),len(z)); k3=float(n3)
            while True:
                na=[max(1,min(nn[i],max(2,int(round(k3*r))))) for i,r in enumerate((1.,0.75,0.625))]
                if na[0]*na[1]*na[2]<=ARROWS_3D_MAX or k3<=2: break
                k3*=0.93
            d=Lm/max(k3,1.)
            g=[np.linspace(c[0]+(c[-1]-c[0])/(2.*k),c[-1]-(c[-1]-c[0])/(2.*k),k) if k>1 else np.array([0.5*(c[0]+c[-1])])
               for c,k in zip((x,y,z),na)]
            def samp(F):   # F (Nz,Ny,Nx) -> (ny,nx,nz) like meshgrid(..., indexing="xy")
                return _lin_axis(_lin_axis(_lin_axis(F,x,g[0],2),y,g[1],1),z,g[2],0).transpose(1,2,0)
            Xa,Ya,Za=np.meshgrid(g[0],g[1],g[2],indexing="xy")
            Ua,Va,Wa=samp(s["Jx"]),samp(s["Jy"]),samp(s["Jz"])
            Ma=np.nan_to_num(np.sqrt(Ua**2+Va**2+Wa**2))
            nar=int(np.count_nonzero(Ma>1e-30))
            if Ma.max()>1e-30:
                msf=np.where(Ma>1e-30,Ma,1.); Lr=0.64*d
                from matplotlib.colors import LogNorm
                mf=Ma.flatten(); mpos=mf[mf>1e-30]
                norm=LogNorm(vmin=max(mpos.min(),Ma.max()*1e-4),vmax=Ma.max())
                cols=plt.cm.viridis(norm(np.clip(mf,norm.vmin,None)))
                lw=1.2*min(1.,(8./max(na))**0.5)
                aQ3.quiver(Xa.flatten(),Za.flatten(),-Ya.flatten(),
                           np.nan_to_num(Ua/msf*Lr).flatten(),np.nan_to_num(Wa/msf*Lr).flatten(),
                           np.nan_to_num(-Va/msf*Lr).flatten(),colors=cols,arrow_length_ratio=0.4,linewidth=lw)
            aQ3.set_xlim(0,x.max()); aQ3.set_ylim(0,z.max()); aQ3.set_zlim(-y.max(),0)
            self._aspect3(aQ3,x.max(),y.max(),z.max())
            sax3(aQ3,"x","z","depth y",f"3-D current field ({nar} arrows)"); _depth_axis(aQ3)
        else:
            x,y,z=s["x"],s["y"],s["z"]
            if not self.cur_st3.get():
                aQ3.text2D(0.5,0.5,"3-D arrows hidden",transform=aQ3.transAxes,ha="center",fontsize=10,color=C["sub"])
            aQ3.set_xlim(0,x.max()); aQ3.set_ylim(0,z.max()); aQ3.set_zlim(-y.max(),0)
            self._aspect3(aQ3,x.max(),y.max(),z.max())
            sax3(aQ3,"x","z","depth y","3-D current field (arrows off)"); _depth_axis(aQ3)
        try: fig.tight_layout(pad=1.)
        except Exception: pass
        self._flush_labels(fig)
        self.cvc.draw_idle()

    def _struct_plane(self,ax,hx,vy,own2,regs,alpha):
        """The structure in a solution plane: each node in its region's colour,
        see-through (alpha, or the part's own opacity from Layers…; hidden
        parts left out) - drawn over the |J| map or alone, under the arrows."""
        if own2 is None or not regs: return
        from matplotlib.colors import ListedColormap,BoundaryNorm,to_rgba
        n=len(regs); tab=np.zeros((n+1,4))
        for k,r in enumerate(regs):
            hid,own=self._lay_part("r",k,r)
            if hid: continue
            c=to_rgba(self._reg_color(r)[0]); a=own if own is not None else alpha
            tab[k]=(c[0],c[1],c[2],max(0.,min(1.,a)))
        idx=np.where(np.asarray(own2)>=0,own2,n)
        ax.pcolormesh(hx,vy,idx,cmap=ListedColormap(tab),norm=BoundaryNorm(np.arange(-0.5,n+1.5,1.),n+1),
                      shading="nearest")
    def _plot_prof(self,ax,s,key):
        """'Current under a contact' panel: J through the chosen contact(s)."""
        dev=s.get("dev") or {}; cons=dev.get("cons") or self.cons
        idx=group_of(cons,key) if key else []
        ttl=f"Current through '{key}'"
        pr=contact_profile(s,cons,idx) if idx else None
        if pr is None:
            msg=(f"'{key}' is not a contact of this solution" if not idx else
                 "Box contacts and gates have no face profile -\nuse the line cut instead")
            ax.text(0.5,0.5,msg,transform=ax.transAxes,ha="center",va="center",fontsize=9,color=C["sub"])
            sax(ax,"","",ttl); return
        st=profile_stats(pr); U,Vn=pr["axes"]; parts=pr["parts"]
        sg=1. if st["mean"]>=0 else -1.
        lbl="J into the device (A/cm²)" if sg>0 else "J out of the device (A/cm²)"
        cu,cv=s[U.lower()],s[Vn.lower()]
        nu=len(np.unique(np.concatenate([p[1] for p in parts]))); nv=max(len(p[2]) for p in parts)
        if nu>=4 and nv>=4 and len(cv)>=6:                       # a contact area: map it
            lo=min(float((sg*p[3]).min()) for p in parts); hi=max(float((sg*p[3]).max()) for p in parts)
            if hi<=lo: hi=lo+1e-30
            pc=None
            for lab,u,v,J,w in parts:
                pc=ax.pcolormesh(u,v,(sg*J).T,cmap="viridis",vmin=lo,vmax=hi,shading="nearest")
            ax.get_figure().colorbar(pc,ax=ax,label=lbl,shrink=0.85,pad=0.02)
            ax.set_xlim(cu[0],cu[-1]); ax.set_ylim(cv[0],cv[-1])
            if Vn=="Y": ax.invert_yaxis()
            if 1/3<=(cu[-1]-cu[0])/max(cv[-1]-cv[0],1e-12)<=3: ax.set_aspect("equal",adjustable="box")
            sax(ax,f"{U.lower()} (nm)",f"{Vn.lower()} (nm)",ttl+" - map over the contact")
        else:                                                    # across the contact
            col=C["green"]
            for n,(lab,u,v,J,w) in enumerate(parts):
                ww=w.sum(axis=1); Ju=sg*(J*w).sum(axis=1)/np.where(ww>0,ww,1.)
                ax.plot(u,Ju,color=col,lw=1.8,marker="o",ms=3,label=lab if len(parts)<=4 else (key if n==0 else None))
            ax.axhline(sg*st["mean"],color=C["yellow"],lw=1,ls="--",label=f"mean {sg*st['mean']:.3g} A/cm²")
            ax.axhline(0,color=C["border"],lw=0.5)
            ax.set_xlim(cu[0],cu[-1]); ax.legend(fontsize=7,framealpha=0.85)
            sax(ax,f"{'depth y' if U=='Y' else U.lower()} (nm)",lbl,ttl)
        txt=f"J max / mean = {st['max']/st['mean'] if sg>0 else st['min']/st['mean']:.2f}"
        if st["ce"]:
            c0,e0=st["ce"][0]
            if e0: txt+=f"\ncentre / edge = {c0/e0:.2f}"
        if pr["skipped"]: txt+=f"\n(not shown: {', '.join(pr['skipped'][:3])})"
        ax.text(0.02,0.97,txt,transform=ax.transAxes,va="top",fontsize=8,color=C["yellow"],
                bbox=dict(boxstyle="round,pad=0.35",facecolor=C["panel"],edgecolor=C["border"],alpha=0.9))

    # ─── FIELD & MOBILITY ──────────────────────────────────────
    def _plot_fi(self):
        s=self.sol; fig=self.ff; fig.clear()
        if s is None: self.cvf.draw_idle(); return
        vis=self._vis("tf")
        if not vis: self._no_panels(fig,self.cvf); return
        axs=row_axes(fig,[(k,1.,None) for k in vis],wspace=0.38)
        aE=axs.get("E") or Figure().add_subplot(); aM=axs.get("mu") or Figure().add_subplot()   # hidden: off screen
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
        self._hmaps={k:v for k,v in getattr(self,"_hmaps",{}).items() if k.figure is not self.fsl}
        self.sl_hover.config(text="")
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
            if qty=="trap density":
                Nt=np.nan_to_num(s.get("Ntrap",np.zeros_like(s["phi"])))
                return np.where(semi&(Nt>0),np.log10(np.maximum(Nt,1e-30)),np.nan),"magma","log₁₀ N_trap, all deep levels (cm⁻³)"
            if qty=="trapped charge":
                Qt=np.nan_to_num(s.get("Qtrap",np.zeros_like(s["phi"])))
                return (np.where(semi,np.sign(Qt)*np.log10(1.+np.abs(Qt)),np.nan),"RdBu_r",
                        "trapped charge ±log₁₀(1+|ρ/q|) (cm⁻³, + = positive)")
            if qty=="log10 G avalanche":
                G=np.nan_to_num(s.get("Gii",np.zeros_like(s["phi"])))
                return (np.where(semi&(G>0),np.log10(np.maximum(G,1e-300)),np.nan),"inferno",
                        "log₁₀ G impact ionisation (cm⁻³s⁻¹)"+("" if self.sol.get("models",0)&M_II else " - avalanche off"))
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
        if qty in ("log10|J gate|","log10 G avalanche"): vmin=max(vmin,vmax-8.)          # the top 8 decades
        if vmin==vmax: vmax=vmin+1e-10
        norm=plt.Normalize(vmin,vmax); sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap)
        fig=self.fsl; fig.clear()
        vis=self._vis("ts")
        if not vis: self._no_panels(fig,self.cvsl); return
        flat=[k for k in vis if k!="3d"]; a3=None; A2={}
        if "3d" in vis and flat:
            gs=GridSpec(len(flat),2,figure=fig,width_ratios=[1.35,1.0],hspace=0.55,wspace=0.18)
            a3=fig.add_subplot(gs[:,0],projection="3d")
            for i_,k in enumerate(flat): A2[k]=fig.add_subplot(gs[i_,1])
        elif "3d" in vis: a3=fig.add_subplot(111,projection="3d")
        else:
            gs=GridSpec(1,len(flat),figure=fig,wspace=0.32)
            for i_,k in enumerate(flat): A2[k]=fig.add_subplot(gs[0,i_])
        for a in A2.values(): a.set_facecolor(C["bg"])
        axy,ayz,axz=A2.get("xy"),A2.get("yz"),A2.get("xz")
        if a3 is not None:
            a3.set_facecolor(C["bg"])
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
        else: fig.colorbar(sm,ax=list(A2.values()),shrink=0.8,pad=0.03,label=lbl)
        own=self._owner(s); dev=s.get("dev") or {}; regs=dev.get("regs",[])
        def pmap(ax,hc,vc,Z,title,xl,yl,hl=None,vl=None,inv=False,o2=None,ovax=None,ovpos=None):
            cuts=[("v",hl),("h",vl)]
            tri_sec=Qr is not None and ovax is not None and ovax[2]=="xyz"[pm.ax]
            if tri_sec:
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
                if inv: ax.invert_yaxis()
                if o2 is not None and dev:
                    hf=np.linspace(hc[0],hc[-1],400); vf=np.linspace(vc[0],vc[-1],400)
                    H2,V2=np.meshgrid(hf,vf)
                    P={ovax[0]:H2,ovax[1]:V2,ovax[2]:np.full_like(H2,ovpos)}
                    o2f=SM.region_index(regs,P["x"],P["y"],P["z"],dev.get("L"))
                    self._ov_plane(ax,hf,vf,o2f,regs,dev.get("cons"),ovax,ovpos,dev.get("L"),cuts)
            else:
                pc=ax.pcolormesh(hc,vc,np.ma.masked_invalid(Z),cmap=cmap,vmin=vmin,vmax=vmax,shading="nearest")
                ax.set_xlim(hc[0],hc[-1]); ax.set_ylim(vc[0],vc[-1])
                if inv: ax.invert_yaxis()
                if o2 is not None:
                    self._ov_plane(ax,hc,vc,o2,regs,dev.get("cons"),ovax,ovpos,dev.get("L"),cuts)
            if hl is not None: ax.axvline(hl,color=C["yellow"],lw=0.9,ls="--",alpha=0.8)
            if vl is not None: ax.axhline(vl,color=C["green"],lw=0.9,ls="--",alpha=0.8)
            ax.set_title(title,fontsize=8,color=C["text"])
            ax.set_xlabel(xl,fontsize=7,color=C["sub"]); ax.set_ylabel(yl,fontsize=7,color=C["sub"])
            ax.tick_params(labelsize=6,colors=C["sub"])
            if ovax is not None: self._reg_hover(ax,hc,vc,Z,lbl,ovax,ovpos)
        if axy is not None:
            pmap(axy,x,y,Zxy,f"XY @ z = {z[sz]:.1f} nm","x (nm)","depth y (nm)",x[sx],y[sy],True,
                 own[sz,:,:] if own is not None else None,("x","y","z"),z[sz])
        if ayz is not None:
            pmap(ayz,z,y,Zyz.T,f"YZ @ x = {x[sx]:.1f} nm","z (nm)","depth y (nm)",z[sz],y[sy],True,
                 own[:,:,sx].T if own is not None else None,("z","y","x"),x[sx])
        if axz is not None:
            pmap(axz,x,z,Zxz,f"XZ @ y = {y[sy]:.1f} nm","x (nm)","z (nm)",x[sx],z[sz],False,
                 own[:,sy,:] if own is not None else None,("x","z","y"),y[sy])
        self._flush_labels(fig)
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
                   Vf=np.array([np.nan if r["Vf"] is None else r["Vf"] for r in recs],float)[o],
                   Qn=np.array([r.get("Qn",np.nan) for r in recs],float)[o],
                   Qp=np.array([r.get("Qp",np.nan) for r in recs],float)[o])
            return dict(R=R,A=None,cons=m["cons"],groups=m["groups"],tkey=m["tkey"],fkey=m["fkey"],
                        okey=m["okey"],bv=m["bv"],wl=m.get("wl"),models=m.get("models"),tname=self.tname)
        return self.ivd

    def _iv_prep(self,d):
        """Terminal (net) currents of a sweep record (the last sweep, a live
        one, or a pinned one) and which terminal to show."""
        R=d["R"]; V=R["V"]; I=R["I"]; F=R["F"]; cv=R["conv"]
        groups=d["groups"]
        bcs=d.get("bcs") or [c.bc for c in d["cons"]]
        gI={k:I[:,idx].sum(axis=1) for k,idx in groups}
        gF={k:F[:,idx].sum(axis=1) for k,idx in groups}
        gate={k:all(bcs[c]==2 for c in idx) for k,idx in groups}
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
        # a floating terminal makes a transfer curve (inverter output) - unless
        # the swept terminal's own current is observed (an open-base transistor)
        return dict(V=V,cv=cv,gI=gI,gF=gF,gA=gA,gate=gate,tkey=tkey,fkey=fkey,okey=okey,cand=cand,
                    groups=groups,R=R,models=d.get("models"),bv=d.get("bv",0.),wl=d.get("wl"),
                    tname=d.get("tname",""),vtc=bool(fkey) and okey!=tkey,T=d.get("T"),vapp=d.get("vapp"))

    @staticmethod
    def _bjt_roles(P):
        """(base, collector, emitter) when the sweep is a Gummel plot - the
        swept terminal drives the base, the observed one is the collector,
        the emitter is the other carrier terminal with the largest current;
        None for anything else (a gain far above 1e5 is not a BJT either)."""
        t,o=P["tkey"],P["okey"]
        if P["fkey"] or t==o or P["gate"].get(t) or P["gate"].get(o) or t not in P["gI"]: return None
        rest=[k for k in P["cand"] if k not in (t,o)]
        if not rest: return None
        gI=P["gI"]; e=max(rest,key=lambda k:float(np.nanmax(np.abs(gI[k]))) if len(gI[k]) else 0.)
        ok=P["cv"]&(np.abs(gI[t])>0)
        if ok.sum()<3: return None
        b=np.abs(gI[o][ok])/np.abs(gI[t][ok])
        if not (0.3<float(np.nanmedian(b))<1e5): return None
        return t,o,e

    def _bjt_curves(self,P):
        """beta and f_T along a Gummel sweep record (see bjt_gain_ft); the
        stored charge is the base's majority carrier: holes when the
        collector current flows in (NPN), electrons otherwise (PNP)."""
        roles=self._bjt_roles(P)
        if roles is None: return None
        b,c,e=roles; R=P["R"]; Ic=P["gI"][c]
        j=int(np.nanargmax(np.abs(Ic))) if len(Ic) else 0
        npn=Ic[j]>=0
        Q=R.get("Qp" if npn else "Qn")
        g=bjt_gain_ft(P["V"],P["gI"][b],Ic,Q,P["cv"])
        if g is None: return None
        g.update(base=b,coll=c,emit=e,npn=npn,AE=P["gA"].get(e,0.))
        return g

    def _plot_iv(self,live=False):
        if not live: self._ptab_fill()
        d=self._iv_data(live); fig=self.fiv
        pins=[p for p in self.pins if p.get("show",True)]
        have=d is not None and len(d["R"]["V"])>0
        if not have and not pins:
            if not live: self._placeholder(fig,self.cviv,"No sweep yet - set it up in the Sweep page and press ↔ Sweep (F6)")
            return
        P=self._iv_prep(d) if have else None
        if P is not None:                       # only pins of the same swept terminal overlay
            pins=[p for p in pins if p["d"]["tkey"]==P["tkey"]]
        self._pins_drawn=len(pins)
        par=[]
        if P is not None and P["vtc"]:
            fig.clear()
            self._plot_vtc(fig,P,par)
        else:
            vis=self._vis("ti")
            if not vis: self._no_panels(fig,self.cviv); return
            fig.clear()
            self._plot_ivc(fig,P,[p for p in pins if not self._iv_prep(p["d"])["vtc"]],par,vis)
        try: fig.tight_layout(pad=1.)
        except Exception: pass
        self.cviv.draw_idle()
        if not live:
            if P is not None:
                hdr=f"{P['tname'] or self.tname}: sweep of {P['tkey']}"+(f", {P['fkey']} floating" if P["fkey"] else "")
            else: hdr="Pinned sweeps"
            self._settext(self.res_par,hdr+"\n"+("\n".join(par) if par else "(no parameters apply to this sweep)"))
    def _replot_iv(self):
        if self.ivd is not None or self.pins: self._plot_iv()

    def _plot_ivc(self,fig,P,pins,par,vis):
        """I-V panels: linear, semi-log (all terminals), BJT gain and f_T;
        P = the last sweep (None: pinned curves only), pins = pinned sweeps
        drawn over it (dashed, in their own colours)."""
        PP=[(p,self._iv_prep(p["d"])) for p in pins]
        bj=self._bjt_curves(P) if P is not None else None
        bjp=[(p,self._bjt_curves(Q)) for p,Q in PP]
        if "bjt" in vis and bj is None and not any(g is not None for _,g in bjp): vis=[k for k in vis if k!="bjt"]
        if not vis: vis=["lin"]
        axs=row_axes(fig,[(k,1.,None) for k in vis],wspace=0.45)
        for a in axs.values(): a.set_facecolor(C["bg"])
        aL=axs.get("lin") or Figure().add_subplot(); aS=axs.get("log") or Figure().add_subplot()
        aB=axs.get("bjt")
        # common current unit for the linear panel
        mags=[]
        if P is not None and len(P["gI"][P["okey"]]): mags.append(float(np.nanmax(np.abs(P["gI"][P["okey"]]))))
        for _,Q in PP:
            if len(Q["gI"][Q["okey"]]): mags.append(float(np.nanmax(np.abs(Q["gI"][Q["okey"]]))))
        sc,un=_cur_unit(max([m for m in mags if np.isfinite(m)] or [0.]))
        aL.axhline(0,color=C["border"],lw=0.5); aL.axvline(0,color=C["border"],lw=0.5)
        aL._iv_sc=sc
        if P is not None:
            self._iv_main(aL,aS,P,sc,un,par,bj)
        for p,Q in PP:
            Io=Q["gI"][Q["okey"]]; ok=Q["cv"]&np.isfinite(Io)
            if not ok.any(): continue
            lab=f"{p['name']}: I({Q['okey']})"
            aL.plot(Q["V"][ok],Io[ok]/sc,color=p["color"],lw=1.6,ls="--",marker="o",ms=2.5,label=lab)
            a=np.abs(Io); okl=ok&(a>0)&(a>=Q["gF"][Q["okey"]])
            if okl.any(): aS.semilogy(Q["V"][okl],a[okl],color=p["color"],lw=1.6,ls="--",marker="o",ms=2.5,label=lab)
        if P is None and PP:
            Q=PP[0][1]
            sax(aL,f"V({Q['tkey']}) (V)",f"I ({un} per cell)","Pinned sweeps")
            sax(aS,f"V({Q['tkey']}) (V)","|I| (A per cell)","Pinned sweeps (semi-log)")
        if aL.get_legend_handles_labels()[0]: aL.legend(fontsize=7)
        if aS.get_legend_handles_labels()[0]: aS.legend(fontsize=6)
        Aok=getattr(aL,"_iv_Aok",0.)
        if Aok>0 and "lin" in axs:
            a2=aL.twinx(); y0,y1=aL.get_ylim(); a2.set_ylim(y0*sc/(Aok*1e4),y1*sc/(Aok*1e4))
            from matplotlib.ticker import EngFormatter
            a2.yaxis.set_major_formatter(EngFormatter(places=None,sep=" "))
            a2.set_ylabel("J (A/cm²)",fontsize=8,color=C["sub"]); a2.tick_params(labelsize=7,colors=C["sub"])
        if aB is not None: self._plot_bjt_panel(aB,P,bj,bjp)
        for p,Q in PP:
            tn=Q["tname"] or ""
            par.append(f"── pinned: {p['name']} - sweep of {Q['tkey']}"+(f" ({tn})" if tn and p["name"] not in tn else "")+" ──")
            self._iv_params(Q,par,None,None,self._bjt_curves(Q))
        n=len(self.pins); sh=len(PP)
        self.iv_pin_lbl.config(text=(f"{n} pinned, {sh} drawn" + ("" if sh==n else
                                     " (hidden, or another swept terminal)")) if n else "")

    def _iv_main(self,aL,aS,P,sc,un,par,bj):
        """The last sweep in detail: the observed terminal on the linear
        panel (converged / below resolution / failed points), every terminal
        on the semi-log one, and the extracted parameters."""
        V,cv,gI,gF,gA,gate=P["V"],P["cv"],P["gI"],P["gF"],P["gA"],P["gate"]
        tkey,okey,cand=P["tkey"],P["okey"],P["cand"]
        Io=gI[okey]; Fo=gF[okey]
        good=cv&(np.abs(Io)>=Fo); low=cv&~(np.abs(Io)>=Fo)
        if good.any(): aL.plot(V[good],Io[good]/sc,color=C["green"],lw=2,marker="o",ms=3.5,label="converged")
        if low.any(): aL.plot(V[low],Io[low]/sc,color=C["sub"],marker="o",mfc="none",ms=4,ls="none",label="below resolution")
        if (~cv).any(): aL.plot(V[~cv],Io[~cv]/sc,color=C["red"],marker="x",ms=6,ls="none",label="not converged")
        sax(aL,f"V({tkey}) (V)",f"I({okey}) ({un} per cell)",f"{okey} current vs {tkey} voltage")
        aL.yaxis.get_major_formatter().set_useOffset(False)
        aL._iv_Aok=gA.get(okey,0.)             # the J axis is added once every curve is drawn
        tun=[k for k in gI if gate.get(k) and len(gI[k]) and np.nanmax(np.abs(gI[k]))>0]
        # the semi-log axis spans the plotted output current and its resolution
        # (+3 decades): a 1e-38 A tunnelling current under a uA drain current
        # would otherwise squeeze the whole subthreshold region into the top
        ao=np.abs(Io); ref=list(ao[good&(ao>0)])+list(np.asarray(gF[okey])[np.asarray(gF[okey])>0])
        ylo=(min(ref)*1e-3) if ref else None
        for n,k in enumerate(cand+[k for k in tun if k not in cand]):
            Ic=np.abs(gI[k]); ok=cv&(Ic>=gF[k])&(Ic>0)
            if ok.any():
                below=ylo is not None and float(Ic[ok].max())<ylo
                aS.semilogy(V[ok],Ic[ok],color=TCOLS[n%len(TCOLS)],lw=2.2 if k==okey else 1.2,
                            ls="--" if k in tun else "-",marker="o" if k==okey else None,ms=3,
                            label=f"|I| {k}"+(" (tunnelling)" if k in tun else "")+
                                  (f" < {ylo:.0e} A, below the axis" if below else ""))
        if np.any(gF[okey]>0):
            aS.semilogy(V,np.maximum(gF[okey],1e-300),color=C["sub"],lw=0.8,ls=":",label=f"resolution ({okey})")
        if ylo is not None and aS.get_ylim()[0]<ylo: aS.set_ylim(bottom=ylo)
        sax(aS,f"V({tkey}) (V)","|I| (A per cell)","All terminals (semi-log)")
        fk=P["fkey"]
        if fk:                                  # an open (floating, I = 0) terminal: its voltage
            Vf=np.asarray(P["R"]["Vf"],float); okf=cv&np.isfinite(Vf)
            if okf.any():
                a3=aS.twinx(); a3.plot(V[okf],Vf[okf],color=C["orange"],lw=1.2,ls="-.",label=f"V({fk}), open")
                a3.set_ylabel(f"V({fk}) - open, I = 0 (V)",fontsize=8,color=C["orange"])
                a3.tick_params(labelsize=7,colors=C["orange"])
        self._iv_params(P,par,aL,aS,bj)

    def _iv_params(self,P,par,aL,aS,bj,vals=None):
        """Parameters extracted from one sweep record: text lines into par,
        the numbers into vals (name -> (value, unit, how to compare: 'd'
        difference, 'p' difference and per cent, 'r' ratio); the parameter
        table) and, with the axes (the last sweep), their markers."""
        V,cv,gI,gF,gA,gate=P["V"],P["cv"],P["gI"],P["gF"],P["gA"],P["gate"]
        tkey,okey,bv,wl=P["tkey"],P["okey"],P["bv"],P["wl"]
        draw=aL is not None
        def put(k,v,u="",c="d",note=""):
            if vals is not None: vals[k]=(v,u,c,note)
        Io=gI[okey]; Fo=gF[okey]; good=cv&(np.abs(Io)>=Fo)
        a=np.abs(Io); okp=good&(a>0)
        box=lambda ax,txt: ax.text(0.03,0.97,txt,transform=ax.transAxes,va="top",fontsize=8,color=C["yellow"],
                                   bbox=dict(boxstyle="round,pad=0.4",facecolor=C["panel"],edgecolor=C["border"],alpha=0.9))
        if gate.get(okey,False) and okp.any():
            j=int(np.argmax(a[okp])); Aa=gA.get(okey,0.)
            Jg=a[okp][j]/(Aa*1e4) if Aa>0 else float("nan")
            par.append(f"Gate tunnelling current |I_G| = {a[okp][j]:.3e} A per cell at V = {V[okp][j]:.3f} V"
                       +(f"  (J_G = {Jg:.3e} A/cm²)" if Jg==Jg else ""))
            put("Gate current |I_G|,max (per cell)",float(a[okp][j]),"A","r")
            if Jg==Jg: put("Gate current density J_G,max",float(Jg),"A/cm²","r")
            if tkey==okey and okp.sum()>=3 and Aa>0:
                o=np.argsort(V[okp]); Vk=V[okp][o]; Jk=a[okp][o]/(Aa*1e4)
                for Vq in (0.5,1.0,-0.5,-1.0):
                    if Vk.min()<=Vq<=Vk.max():
                        Jq=float(10**np.interp(Vq,Vk,np.log10(Jk)))
                        par.append(f"J_G({Vq:+.1f} V) = {Jq:.3e} A/cm²"); put(f"J_G at {Vq:+.1f} V",Jq,"A/cm²","r")
        if gate.get(tkey,False) and not gate.get(okey,False) and okp.sum()>=3:
            n0=len(par)
            Vg=V[okp]; Ia=a[okp]
            sw=subthreshold_swing(Vg,Ia,floor=Fo[okp])
            if sw is not None:
                SS=sw["SS"]; par.append(f"Subthreshold swing SS = {SS:.1f} mV/dec (steepest decade, "
                                        f"{sw['V0']:.3g} … {sw['V1']:.3g} V)")
                put("Subthreshold swing SS",SS,"mV/dec","d")
            # drain-source voltage: observed terminal against the farthest other carrier terminal
            Vds=None; vapp=P.get("vapp")
            if vapp is not None:
                gd=dict(P["groups"])
                vk=lambda k: float(np.mean([vapp[c] for c in gd[k]]))
                oth=[k for k in gd if not gate.get(k,False) and k not in (okey,tkey)]
                if oth:
                    vo=vk(okey); vs=[vk(k) for k in oth]
                    Vds = vo-min(vs) if vo>=min(vs) else vo-max(vs)
            th=threshold_voltage(Vg,Ia,Vds)
            if th is not None:
                Vth=th["Vth"]; edge=th["edge"]
                how=("linear extrapolation at max g_m"+(f", V_DS/2 = {0.5*Vds*1e3:.0f} mV subtracted" if Vds else "")
                     if th["method"]=="lin" else f"saturation, √I_D extrapolation at V_DS = {Vds:.3g} V")
                par.append(f"V_th ({how}) = {Vth:.3f} V"+
                           ("  - still steepening at the sweep end: extend the sweep for a firm value" if edge else ""))
                par.append(f"g_m,max = {th['gm']:.3e} S per cell at V = {th['Vgm']:.3f} V")
                put("V_th (extrapolation)" if th["method"]=="lin" else "V_th (√I_D extrapolation, saturation)",
                    Vth,"V","d","still steepening at the sweep end" if edge else "")
                put("g_m,max (per cell)",th["gm"],"S","r")
                if draw:
                    aL.axvline(Vth,color=C["yellow"],lw=1,ls="--",alpha=0.7)
                    aL.text(Vth,0.98,f" V_th {Vth:.2f} V",color=C["yellow"],fontsize=8,va="top",
                            transform=aL.get_xaxis_transform())
            if wl:
                Icc=1e-7*wl; k=np.where((Ia[:-1]<Icc)&(Ia[1:]>=Icc)|(Ia[:-1]>=Icc)&(Ia[1:]<Icc))[0]
                if len(k):
                    k=k[0]; l0,l1=np.log(Ia[k]),np.log(Ia[k+1])
                    Vcc=Vg[k]+(Vg[k+1]-Vg[k])*(np.log(Icc)-l0)/(l1-l0)
                    par.append(f"V_th (constant current 100 nA x W/L = {Icc:.2e} A) = {Vcc:.3f} V")
                    put("V_th (constant current, 100 nA·W/L)",float(Vcc),"V","d")
            Ion,Ioff=Ia[-1],Ia[0]
            if abs(Vg[0])>abs(Vg[-1]): Ion,Ioff=Ioff,Ion
            par.append(f"I_on = {Ion:.3e} A, I_off = {Ioff:.3e} A per cell"+(f"  (I_on/I_off = {Ion/Ioff:.2e})" if Ioff>0 else ""))
            put("I_on (per cell)",float(Ion),"A","r"); put("I_off (per cell)",float(Ioff),"A","r")
            if Ioff>0: put("I_on / I_off",float(Ion/Ioff),"","r")
            if draw and len(par)>n0: box(aS,"\n".join(p.split("  - ")[0] for p in par[n0:n0+2]))
        Aok_cm2=gA.get(okey,0.)*1e4
        if tkey==okey and Aok_cm2>0 and not gate.get(okey,False):
            J=Io/Aok_cm2
            fwd=good&(J>1e-20)
            dio=diode_ideality(V[fwd],J[fwd],T=P.get("T") or 300.,floor=(Fo/Aok_cm2)[fwd])
            if dio is not None and 0.8<dio["n"]<4.:
                n_id,J0=dio["n"],dio["J0"]
                par.append(f"Diode ideality n = {n_id:.3f}, J0 = {J0:.2e} A/cm² (straight part {dio['V0']:.3g} … "
                           f"{dio['V1']:.3g} V, kT/q at {P.get('T') or 300.:g} K)")
                put("Diode ideality factor n",float(n_id),"","d"); put("Saturation current density J_0",float(J0),"A/cm²","r")
                if draw: box(aS,f"n = {n_id:.3f}\nJ₀ = {J0:.2e} A/cm²")
        # output characteristic of a transistor: Early voltage
        if tkey==okey and not gate.get(okey,False) and len(P["groups"])>=3:
            ev=early_voltage(V,Io,cv)
            if ev is not None and ev["VA"]<1e5:
                par.append(f"Early voltage V_A = {ev['VA']:.3g} V (output conductance g_o = {ev['go']:.3e} S per cell; "
                           f"straight line over {ev['V0']:.3g} … {ev['V1']:.3g} V)")
                put("Early voltage V_A",float(ev["VA"]),"V","p"); put("Output conductance g_o (per cell)",float(ev["go"]),"S","r")
                if draw:
                    xv=np.linspace(V.min(),V.max(),50); sc=getattr(aL,"_iv_sc",None) or 1.
                    aL.plot(xv,(ev["a"]+ev["b"]*xv)/sc,color=C["yellow"],lw=1,ls=":",
                            label=f"V_A = {ev['VA']:.3g} V (line meets I = 0 at -V_A)")
        # breakdown
        mdl=P.get("models")
        aval=mdl is not None and bool(int(mdl)&M_II)
        if Aok_cm2>0 and not gate.get(okey,False) and np.max(np.abs(V))>=3.:
            Ja=np.abs(Io)/Aok_cm2                    # (failed points: their partial state)
            Vbd=None; how=""
            if aval:
                Vbd,how=avalanche_bv(V,Ja,cv,thr=bv,floor=np.asarray(Fo)/Aok_cm2)
                if Vbd is not None:
                    par.append(f"Avalanche breakdown BV = {abs(Vbd):.4g} V ({how}; impact ionisation on)")
                    put("Breakdown voltage BV (avalanche)",float(abs(Vbd)),"V","p")
                else:
                    Vm=float(V[np.argmax(np.abs(V))])
                    par.append(f"No avalanche breakdown up to V = {Vm:.4g} V (impact ionisation on)")
                    put("Breakdown voltage BV (avalanche)",f"> {abs(Vm):.4g} V","V","")
            elif bv>0:
                for k in range(1,len(Ja)):
                    if good[k] and good[k-1] and Ja[k-1]>1e-30 and Ja[k]/Ja[k-1]>10. and Ja[k]>bv: Vbd=V[k]; break
                if Vbd is None:
                    for k in range(len(Ja)):
                        if good[k] and Ja[k]>bv: Vbd=V[k]; break
                if Vbd is not None:
                    par.append(f"Current knee at V = {Vbd:.2f} V (|J| > {bv:g} A/cm²; heuristic - impact ionisation "
                               f"is off: Solver page)")
                    put(f"Current knee (|J| > {bv:g} A/cm², avalanche off)",float(Vbd),"V","p")
            if draw and Vbd is not None:
                for ax in (aL,aS): ax.axvline(Vbd,color="#ff4444",lw=2,ls="--",alpha=0.8)
                aL.text(Vbd,0.98,f" BV {abs(Vbd):.4g} V" if aval else f" knee {Vbd:.3g} V",
                        color="#ff6666",fontsize=8,va="top",transform=aL.get_xaxis_transform())
        # Gummel plot of a BJT: gain and transit frequency
        if bj is not None:
            Ic=bj["Ic"]; beta=bj["beta"]; k=int(np.nanargmax(beta))
            par.append(f"Current gain β = I_C/I_B: {beta[k]:.1f} at most (I_C = {Ic[k]:.3e} A per cell), "
                       f"{beta[-1]:.1f} at I_C = {Ic[-1]:.3e} A")
            put("Current gain β,max",float(beta[k]),"","p"); put("I_C at β,max (per cell)",float(Ic[k]),"A","r")
            hk=np.where((np.arange(len(beta))>k)&(beta<0.5*beta[k]))[0]
            if len(hk):
                par.append(f"β halves (high injection, Kirk effect, series resistance) by I_C = {Ic[hk[0]]:.3e} A")
                put("I_C where β has halved (per cell)",float(Ic[hk[0]]),"A","r")
            fT=bj["fT"]; IcT=bj["IcT"]
            if np.isfinite(fT).any():
                j=int(np.nanargmax(fT)); AE=bj["AE"]
                Jc=f", J_C = {IcT[j]/(AE*1e4):.3e} A/cm² over the {bj['emit']} contact" if AE>0 else ""
                par.append(f"f_T,max = {fT[j]/1e9:.3g} GHz at I_C = {IcT[j]:.3e} A{Jc} (quasi-static ΔI_C / 2π ΔQ)")
                put("Transit frequency f_T,max",float(fT[j]),"Hz","p","still rising at the sweep end" if j==len(fT)-1 else "")
                put("I_C at f_T,max (per cell)",float(IcT[j]),"A","r")
                if AE>0: put("J_C at f_T,max",float(IcT[j]/(AE*1e4)),"A/cm²","r")
                if j<len(fT)-1 and np.nanmin(fT[j:])<0.7*fT[j]:
                    par.append("   past the f_T peak the base widens into the collector (Kirk effect) and "
                               "high injection sets in")
                elif j==len(fT)-1:
                    par.append("   f_T still rising at the sweep end - extend the sweep to see its peak")

    def _plot_bjt_panel(self,ax,P,bj,bjp):
        """β (left axis) and f_T (right axis) against I_C, log scale."""
        a2=ax.twinx(); anyft=False; Ic_all=[]
        def draw(g,col,ls,lab,lw):
            nonlocal anyft
            ok=g["Ic"]>0
            ax.semilogx(g["Ic"][ok],g["beta"][ok],color=col,lw=lw,ls=ls,marker="o",ms=2.5,label=f"β {lab}")
            Ic_all.extend(list(g["Ic"][ok]))
            f=g["fT"]; IT=g["IcT"]; okf=np.isfinite(f)&(f>0)&(IT>0)
            if okf.any():
                anyft=True
                a2.semilogx(IT[okf],f[okf]/1e9,color=col,lw=lw,ls=":" if ls=="-" else "-.",label=f"f_T {lab}")
                j=int(np.nanargmax(np.where(okf,f,-1.)))
                a2.plot([IT[j]],[f[j]/1e9],marker="*",ms=11,color=col,ls="none")
        if bj is not None: draw(bj,C["green"],"-","(last sweep)",2.)
        for p,g in bjp:
            if g is not None: draw(g,p["color"],"--",p["name"],1.5)
        if bj is not None and np.isfinite(bj["fT"]).any():
            j=int(np.nanargmax(bj["fT"]))
            ax.axvline(bj["IcT"][j],color=C["yellow"],lw=1,ls="--",alpha=0.7)
            ax.text(bj["IcT"][j],0.98,f" f_T,max {bj['fT'][j]/1e9:.3g} GHz",color=C["yellow"],fontsize=8,va="top",
                    transform=ax.get_xaxis_transform())
        sax(ax,"I_C (A per cell)","β = I_C / I_B","BJT: gain β (left axis) and f_T (right axis) vs I_C")
        a2.set_ylabel("f_T (GHz, quasi-static)" if anyft else "f_T: no stored-charge data",fontsize=8,color=C["orange"])
        a2.tick_params(labelsize=7,colors=C["orange"])
        if anyft: a2.set_ylim(bottom=0)
        ax.set_ylim(bottom=0)
        h1,l1=ax.get_legend_handles_labels(); h2,l2=a2.get_legend_handles_labels()
        if h1 or h2: ax.legend(h1+h2,l1+l2,fontsize=6,loc="lower left",framealpha=0.85)

    def _vtc_params(self,P,par,vals=None):
        """Inverter parameters of a sweep with a floating output (a voltage
        transfer curve): V_M, gain, V_OH/V_OL, V_IL/V_IH, noise margins and
        the peak supply current - text into par, numbers into vals."""
        def put(k,v,u="",c="d",note=""):
            if vals is not None: vals[k]=(v,u,c,note)
        V,cv,gI,okey=P["V"],P["cv"],P["gI"],P["okey"]
        Vf=np.asarray(P["R"]["Vf"],float)
        ok=cv&np.isfinite(Vf)
        Vi,idx=np.unique(V[ok],return_index=True); Vo=Vf[ok][idx]
        r=dict(Vi=Vi,Vo=Vo,g=None,VM=None,VIL=None,VIH=None)
        if len(Vi)>=3:
            g=np.gradient(Vo,Vi); r["g"]=g
            d_=Vo-Vi; k=np.where(np.diff(np.sign(d_))!=0)[0]
            if len(k):
                k=k[0]; VM=Vi[k]+(Vi[k+1]-Vi[k])*d_[k]/(d_[k]-d_[k+1]); r["VM"]=VM
                par.append(f"Switching threshold V_M = {VM:.3f} V"); put("Switching threshold V_M",float(VM),"V","d")
            gi=int(np.argmax(np.abs(g))); par.append(f"Max |gain| = {abs(g[gi]):.1f} at V_in = {Vi[gi]:.3f} V")
            put("Max |gain| dV_out/dV_in",float(abs(g[gi])),"","p")
            VOH,VOL=(Vo[0],Vo[-1]) if Vo[0]>Vo[-1] else (Vo[-1],Vo[0])
            par.append(f"V_OH = {VOH:.4f} V, V_OL = {VOL:.4f} V")
            put("Output high V_OH",float(VOH),"V","d"); put("Output low V_OL",float(VOL),"V","d")
            c=np.where(np.diff(np.sign(g+1.))!=0)[0]
            if len(c)>=2:
                f=lambda j:Vi[j]+(Vi[j+1]-Vi[j])*(-1.-g[j])/(g[j+1]-g[j])
                VIL,VIH=f(c[0]),f(c[-1]); r["VIL"],r["VIH"]=VIL,VIH
                par.append(f"V_IL = {VIL:.3f} V, V_IH = {VIH:.3f} V (gain = -1)")
                par.append(f"Noise margins NM_L = {VIL-VOL:.3f} V, NM_H = {VOH-VIH:.3f} V")
                put("Input low V_IL (gain = -1)",float(VIL),"V","d"); put("Input high V_IH (gain = -1)",float(VIH),"V","d")
                put("Noise margin NM_L",float(VIL-VOL),"V","d"); put("Noise margin NM_H",float(VOH-VIH),"V","d")
        Ia=np.abs(gI[okey])[cv]
        if len(Ia):
            j=int(np.argmax(Ia)); par.append(f"Peak |I({okey})| = {Ia[j]:.3e} A per cell at V_in = {V[cv][j]:.3f} V")
            put(f"Peak supply current |I| (per cell)",float(Ia[j]),"A","r")
        return r

    def _plot_vtc(self,fig,P,par):
        V,cv,gI,gF,gate=P["V"],P["cv"],P["gI"],P["gF"],P["gate"]
        tkey,fkey,okey=P["tkey"],P["fkey"],P["okey"]
        Vf=np.asarray(P["R"]["Vf"],float)
        gs=GridSpec(1,2,figure=fig,wspace=0.45)
        aT=fig.add_subplot(gs[0]); aI=fig.add_subplot(gs[1])
        aT.set_facecolor(C["bg"]); aI.set_facecolor(C["bg"])
        ok=cv&np.isfinite(Vf)
        r=self._vtc_params(P,par)
        Vi,Vo,g=r["Vi"],r["Vo"],r["g"]
        aT.plot(V[ok],Vf[ok],color=C["green"],lw=2,marker="o",ms=3,label=f"V({fkey})")
        if (~cv).any(): aT.plot(V[~cv],Vf[~cv],color=C["red"],marker="x",ls="none",ms=6,label="not converged")
        if len(V): aT.plot([V.min(),V.max()],[V.min(),V.max()],color=C["sub"],ls=":",lw=1,label="V_out = V_in")
        if g is not None:
            a2=aT.twinx()
            a2.plot(Vi,g,color=C["orange"],lw=1.2,ls="--",label="gain")
            a2.set_ylabel("gain dV_out/dV_in",fontsize=8,color=C["orange"]); a2.tick_params(labelsize=7,colors=C["orange"])
            if r["VM"] is not None:
                VM=r["VM"]
                aT.axvline(VM,color=C["yellow"],lw=1,ls="--",alpha=0.8)
                aT.text(VM,max(Vo),f" V_M = {VM:.3f} V",color=C["yellow"],fontsize=8,va="top")
            if r["VIL"] is not None:
                for vv in (r["VIL"],r["VIH"]): aT.axvline(vv,color=C["sub"],lw=0.8,ls=":")
        aT.legend(fontsize=7,loc="center left")
        sax(aT,f"V({tkey}) (V)",f"V({fkey}) (V)","Voltage transfer curve (output floating, I = 0)")
        for n,k in enumerate([k for k in gI if not gate[k] and k!=fkey]):
            Ic=np.abs(gI[k]); okk=cv&(Ic>=gF[k])&(Ic>0)
            if okk.any():
                aI.semilogy(V[okk],Ic[okk],color=TCOLS[n%len(TCOLS)],lw=2.2 if k==okey else 1.2,
                            marker="o" if k==okey else None,ms=3,label=f"|I| {k}")
        if np.any(gF[okey]>0): aI.semilogy(V,np.maximum(gF[okey],1e-300),color=C["sub"],lw=0.8,ls=":",label="resolution")
        aI.legend(fontsize=7); sax(aI,f"V({tkey}) (V)","|I| (A per cell)","Supply / terminal currents")

    # ─── PARAMETER TABLE ───────────────────────────────────────
    def _sweep_vals(self,d,P=None):
        """The parameter table's numbers for one sweep record (cached on
        pinned records, which never change)."""
        if P is None: P=self._iv_prep(d)
        vals={}; par=[]
        try:
            if P["vtc"]: self._vtc_params(P,par,vals)
            else: self._iv_params(P,par,None,None,self._bjt_curves(P),vals)
        except Exception:
            self._log("parameter table: "+traceback.format_exc().strip().splitlines()[-1])
        Io=np.asarray(P["gI"][P["okey"]],float); ok=P["cv"]&np.isfinite(Io)
        if ok.any():
            j=int(np.where(ok)[0][-1])
            vals[f"|I_out| at the sweep end (per cell)"]=(float(abs(Io[j])),"A","r","")
        return vals,P

    def _ptab_cols(self):
        """[(name, colour, prepared record, values)] of the last sweep and the
        pinned sweeps the table shows."""
        out=[]
        allp=bool(self.ptab_all.get())
        for p in self.pins:
            if not (allp or p.get("show",True)): continue
            if "_vals" not in p: p["_vals"],p["_P"]=self._sweep_vals(p["d"])
            out.append((p["name"],p["color"],p["_P"],p["_vals"]))
        d=self.ivd
        if d is not None and len(d["R"]["V"]):
            if "_ptab" not in d: d["_ptab"]=self._sweep_vals(d)        # (a new sweep is a new record)
            vals,P=d["_ptab"]
            out.append(("Last sweep",C["green"],P,vals))
        return out

    @staticmethod
    def _ptab_fmt(v,u):
        """A table value with an SI prefix where that reads better."""
        if isinstance(v,str): return v
        if v is None or not np.isfinite(v): return "–"
        if u in ("A","S","Hz"):
            a=abs(v)
            if a==0: return f"0 {u}"
            for f_,pre in ((1e12,"T"),(1e9,"G"),(1e6,"M"),(1e3,"k"),(1.,""),(1e-3,"m"),(1e-6,"µ"),(1e-9,"n"),
                           (1e-12,"p"),(1e-15,"f"),(1e-18,"a")):
                if a>=f_*0.99995: return f"{v/f_:.4g} {pre}{u}"
            return f"{v:.3e} {u}"
        if u=="A/cm²": return f"{v:.3e} {u}"
        if u=="V": return f"{v:.4g} V"
        if u: return f"{v:.4g} {u}"
        return f"{v:.4g}"

    @classmethod
    def _ptab_diff(cls,v,r,u,how):
        """The comparison of v with the reference value r: a difference
        (mV for small voltage differences), with per cent, or a ratio."""
        if isinstance(v,str) or isinstance(r,str) or v is None or r is None: return ""
        if not (np.isfinite(v) and np.isfinite(r)) or not how: return ""
        if how=="r":
            if r==0 or v/r<=0: return ""
            q=v/r; return f"×{q:.3g}" if 1e-3<=q<1e4 else f"×{q:.2e}"
        dv=v-r
        if u=="V" and abs(dv)<1.: txt=f"{dv*1e3:+.4g} mV"
        elif u in ("A","S","Hz"): txt=("+" if dv>=0 else "")+cls._ptab_fmt(dv,u)
        elif u: txt=f"{dv:+.4g} {u}"
        else: txt=f"{dv:+.4g}"
        if how=="p" and r!=0: txt+=f" ({100.*dv/abs(r):+.3g}%)"
        return txt

    def _ptab_rows(self):
        """Header, rows of display text and the raw table (for Copy and CSV)."""
        cols=self._ptab_cols()
        names=[c[0] for c in cols]
        ref=self.ptab_ref.get()
        if ref not in names: ref=names[0] if names else ""
        ri=names.index(ref) if ref in names else -1
        mode=self.ptab_mode.get()
        info=[]
        for nm,col,P,vals in cols:
            V=P["V"]; n=len(V); nc=int(np.sum(P["cv"]))
            sw=f"V({P['tkey']}) {V.min():.3g} … {V.max():.3g} V" if n else P["tkey"]
            if P["fkey"]: sw+=f", {P['fkey']} open"
            T=P.get("T"); mdl=P.get("models")
            if mdl is None: ml="-"
            else:                       # relative to the Solver page's defaults: short, and shows what differs
                mdl=int(mdl); mn=self.MODEL_NAMES
                add=[v for b,v in mn.items() if mdl&b and not M_DEFAULT&b]
                rem=[v for b,v in mn.items() if M_DEFAULT&b and not mdl&b]
                if len(rem)>2: ml=", ".join(v for b,v in mn.items() if mdl&b) or "none"
                else: ml="default"+(" + "+", ".join(add) if add else "")+(" − "+", ".join(rem) if rem else "")
            info.append((P["tname"] or "-",sw,P["okey"],"-" if T is None else f"{T:g} K",ml,f"{nc}/{n}"))
        keys=[]
        for _,_,_,vals in cols:
            for k in vals:
                if k not in keys: keys.append(k)
        rows=[]
        for t,lab in enumerate(("Device","Sweep","Observed terminal","Temperature","Physics models","Converged points")):
            rows.append(("info",lab,"",[i[t] for i in info],[i[t] for i in info]))
        for k in keys:
            ent=[c[3].get(k) for c in cols]
            u=next((e[1] for e in ent if e is not None),"")
            how=next((e[2] for e in ent if e is not None and e[2]),"")
            rv=ent[ri][0] if ri>=0 and ent[ri] is not None else None
            disp=[]; raw=[]
            for j,e in enumerate(ent):
                if e is None: disp.append("–"); raw.append(""); continue
                v,_,_,note=e
                txt=self._ptab_fmt(v,u)+(" *" if note else "")
                raw.append(v)
                if j!=ri and mode!="Values":
                    dtx=self._ptab_diff(v,rv,u,how)
                    if mode=="Differences only": txt=dtx or txt
                    elif dtx: txt=f"{txt}   ({dtx})"
                disp.append(txt)
            rows.append(("par",k,u,disp,raw,ent))
        return cols,names,ri,rows

    def _ptab_win(self):
        top=self._ptab_top
        if top is not None and top.winfo_exists():
            top.deiconify(); top.lift(); self._ptab_fill(); return
        top=tk.Toplevel(self); top.title("Sweep parameters - last sweep and pinned sweeps"); top.configure(bg=C["panel"])
        top.transient(self); self._ptab_top=top
        try: top.geometry(f"1000x520+{self.winfo_rootx()+90}+{self.winfo_rooty()+90}")
        except Exception: pass
        tk.Label(top,text="Parameters extracted from each sweep (as in the I-V tab's text), side by side. Differences are "
                          "taken against the reference column: voltages in mV or V, currents and frequencies as ratios "
                          "(×), gain and breakdown also in per cent. * = see the I-V tab's note (e.g. still rising at the "
                          "sweep end). Pin a sweep (I-V tab) to add a column.",
                 bg=C["panel"],fg=C["sub"],font=self.F["small"],justify="left",anchor="w",wraplength=960
                 ).pack(fill=tk.X,padx=10,pady=(8,2))
        cf=tk.Frame(top,bg=C["panel"]); cf.pack(fill=tk.X,padx=8,pady=2)
        tk.Label(cf,text="Reference",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(2,2))
        self.ptab_rcb=ttk.Combobox(cf,textvariable=self.ptab_ref,values=[],state="readonly",width=28,font=self.F["ui"])
        self.ptab_rcb.pack(side=tk.LEFT,padx=2); self.ptab_rcb.bind("<<ComboboxSelected>>",lambda e:self._ptab_fill())
        tk.Label(cf,text="Show",bg=C["panel"],fg=C["sub"],font=self.F["ui"]).pack(side=tk.LEFT,padx=(10,2))
        cb=ttk.Combobox(cf,textvariable=self.ptab_mode,values=["Values","Values and differences","Differences only"],
                        state="readonly",width=22,font=self.F["ui"])
        cb.pack(side=tk.LEFT,padx=2); cb.bind("<<ComboboxSelected>>",lambda e:self._ptab_fill())
        self._check(cf,"Hidden pins too",self.ptab_all,self._ptab_fill)
        lf=tk.Frame(top,bg=C["panel"]); lf.pack(fill=tk.BOTH,expand=True,padx=8,pady=2)
        tv=ttk.Treeview(lf,columns=("u",),show="tree headings",height=16,selectmode="extended")
        vs=ttk.Scrollbar(lf,orient="vertical",command=tv.yview); hs=ttk.Scrollbar(lf,orient="horizontal",command=tv.xview)
        tv.configure(yscrollcommand=vs.set,xscrollcommand=hs.set)
        tv.grid(row=0,column=0,sticky="nsew"); vs.grid(row=0,column=1,sticky="ns"); hs.grid(row=1,column=0,sticky="ew")
        lf.rowconfigure(0,weight=1); lf.columnconfigure(0,weight=1)
        tv.tag_configure("info",foreground=C["sub"]); tv.tag_configure("odd",background="#0e1a2b")
        self._ptab_tv=tv
        bf=tk.Frame(top,bg=C["panel"]); bf.pack(fill=tk.X,padx=8,pady=(4,8))
        self._btn(bf,"Copy (tab-separated)",self._ptab_copy).pack(side=tk.LEFT,padx=3)
        self._btn(bf,"Export CSV…",self._ptab_csv,bg=C["border"]).pack(side=tk.LEFT,padx=3)
        self._btn(bf,"Pins…",self._pins_win,bg=C["border"]).pack(side=tk.LEFT,padx=3)
        self.ptab_note=tk.Label(bf,text="",bg=C["panel"],fg=C["sub"],font=self.F["small"],anchor="w")
        self.ptab_note.pack(side=tk.LEFT,padx=8,fill=tk.X,expand=True)
        self._btn(bf,"Close",top.destroy,bg=C["border"]).pack(side=tk.RIGHT,padx=3)
        self._ptab_fill()

    def _ptab_fill(self):
        top=getattr(self,"_ptab_top",None)
        if top is None or not top.winfo_exists(): return
        tv=self._ptab_tv
        cols,names,ri,rows=self._ptab_rows()
        self.ptab_rcb.configure(values=names)
        if names and self.ptab_ref.get() not in names: self.ptab_ref.set(names[ri if ri>=0 else 0])
        ids=[f"c{j}" for j in range(len(cols))]
        tv.configure(columns=["u"]+ids)
        tv.heading("#0",text="Parameter",anchor="w"); tv.column("#0",width=310,minwidth=160,stretch=False)
        tv.heading("u",text="Unit"); tv.column("u",width=62,minwidth=40,anchor="center",stretch=False)
        wid=200 if self.ptab_mode.get()=="Values and differences" else 130
        for j,(nm,col,_,_) in enumerate(cols):
            tv.heading(ids[j],text=("★ " if j==ri else "")+nm)
            tv.column(ids[j],width=wid,minwidth=80,anchor="w",stretch=True)
        for it in tv.get_children(): tv.delete(it)
        k=0
        for row in rows:
            kind,lab,u,disp=row[0],row[1],row[2],row[3]
            tags=("info",) if kind=="info" else (("odd",) if k%2 else ())
            if kind!="info": k+=1
            tv.insert("","end",text=lab,values=[u]+disp,tags=tags)
        npar=sum(1 for r in rows if r[0]=="par")
        if not cols:
            self.ptab_note.config(text="No sweep yet: run one (↔ Sweep, F6) and pin it, or pin several to compare.")
        else:
            self.ptab_note.config(text=f"{len(cols)} sweep{'s' if len(cols)!=1 else ''}, {npar} parameters"
                                       +(f" · reference: {names[ri]}" if ri>=0 and len(cols)>1 else ""))

    def _ptab_table(self):
        """The table as rows of raw values (numbers in base units) with the
        comparison against the reference: for Copy and the CSV file."""
        cols,names,ri,rows=self._ptab_rows()
        hdr=["Parameter","Unit"]+names
        others=[j for j in range(len(names)) if j!=ri] if ri>=0 else []
        hdr+=[f"{names[j]} vs {names[ri]}" for j in others]
        out=[hdr]
        for row in rows:
            if row[0]=="info": out.append([row[1],""]+[str(x) for x in row[4]]+[""]*len(others)); continue
            _,lab,u,disp,raw,ent=row
            line=[lab,u]+[("%.6g"%x if isinstance(x,float) else str(x)) for x in raw]
            rv=ent[ri][0] if ri>=0 and ent[ri] is not None else None
            for j in others:
                e=ent[j]
                line.append(self._ptab_diff(e[0],rv,u,e[2]) if e is not None and rv is not None else "")
            out.append(line)
        return out

    def _ptab_copy(self):
        t=self._ptab_table()
        txt="\n".join("\t".join(str(x) for x in r) for r in t)
        self.clipboard_clear(); self.clipboard_append(txt)
        self.ptab_note.config(text=f"Copied {len(t)-1} rows - paste into a spreadsheet")
    def _ptab_csv(self):
        p=filedialog.asksaveasfilename(defaultextension=".csv",filetypes=[("CSV","*.csv")],parent=self._ptab_top,
                                       initialfile="sweep_parameters.csv")
        if not p: return
        import csv
        with open(p,"w",newline="",encoding="utf-8") as f: csv.writer(f).writerows(self._ptab_table())
        self.sv.set(f"Parameter table → {p}"); self._log(f"Parameter table exported to {p}")

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
            d=self._device_json()
            with open(p,"w") as f: json.dump(d,f,indent=1)
            self.sv.set(f"Device saved → {p}"); self._log(f"Device saved to {p}")
        except Exception:
            messagebox.showerror("Save device",traceback.format_exc())
    def _device_json(self):
        """The editor's device as the device-file dict (what Save writes)."""
        return dict(format="semisim3d-device",version=1,name=self.tname,
                   Lx=self.Lx.get(),Ly=self.Ly.get(),Lz=self.Lz.get(),Nx=self.Nx.get(),Ny=self.Ny.get(),Nz=self.Nz.get(),
                   graded=self.use_graded.get(),dense_x=self.dense_x,dense_y=self.dense_y,dense_z=self.dense_z,
                   meta=self.tmeta,regs=[obj_dict(r) for r in self.regs],cons=[obj_dict(c) for c in self.cons],
                   sheets=[obj_dict(s) for s in self.sheets],defects=[d.to_dict() for d in self.defects],
                   mesh=dict(type=self.mesh_type.get(),prism_axis=self.prism_ax.get(),amr=self.amr_nodes or {}),
                   sweep=dict(tgt=self._gmap.get(self.iv_tgt.get()),flt=self._gmap.get(self.iv_flt.get()),
                              out=self._gmap.get(self.iv_out.get()),Vs=self.ivS.get(),Ve=self.ivE.get(),
                              N=self.ivN.get(),bv=self.bv_thr.get()),
                   solver=dict(T=self.T.get(),max_p=self.max_p.get(),max_c=self.max_c.get(),max_g=self.max_g.get(),
                               tol_p=self.tol_p.get(),tol_c=self.tol_c.get(),bgn=self.m_bgn.get(),
                               tau=self.m_tau.get(),vsat=self.m_vsat.get(),qc=self.m_qc.get(),
                               lombardi=self.m_lb.get(),tunnel=self.m_tun.get(),
                               metal_te=self.m_te.get(),metal_tunnel=self.m_sbt.get(),
                               avalanche=self.m_ii.get(),tau0=self.tau0.get()))
    def _open_device(self):
        if self.busy: return
        p=filedialog.askopenfilename(filetypes=[("SemiSim device","*.json"),("All files","*")])
        if not p: return
        try:
            with open(p) as f: d=json.load(f)
            if not isinstance(d,dict) or d.get("format")!="semisim3d-device": raise ValueError("not a SemiSim 3D device file")
            warn=check_device_json(d)
            regs=[Reg(**r) for r in d["regs"]]; cons=[Con(**c) for c in d["cons"]]
            sheets=[Sheet(**s) for s in d.get("sheets",[])]
            defects=[SD.Defect.from_dict(x) for x in d.get("defects",[])]
            bad=[x.get("kind") for x in d.get("defects",[]) if x.get("kind") not in SD.LIB]
            if bad: raise ValueError(f"unknown defect kind(s): {', '.join(map(str,sorted(set(bad))))}")
        except Exception as e:
            messagebox.showerror("Open device",f"Cannot read {p}:\n{e}"); return
        for w in warn: self._log(f"Open device: {w}")
        self._apply_device_json(d,regs,cons,sheets,defects,src=p)
    def _device_parts(self,d):
        """Regions, contacts, sheets and defects of a device-file dict."""
        regs=[Reg(**r) for r in d["regs"]]; cons=[Con(**c) for c in d["cons"]]
        sheets=[Sheet(**s) for s in d.get("sheets",[])]
        defects=[SD.Defect.from_dict(x) for x in d.get("defects",[])]
        return regs,cons,sheets,defects
    def _apply_device_json(self,d,regs,cons,sheets,defects,src="",label=None):
        """Put a device-file dict into the editor (Open device, A/B 'Load')."""
        p=src
        self.tname=d.get("name") or os.path.basename(p); self.tmeta=copy.deepcopy(d.get("meta") or {})
        for v,k in ((self.Lx,"Lx"),(self.Ly,"Ly"),(self.Lz,"Lz"),(self.Nx,"Nx"),(self.Ny,"Ny"),(self.Nz,"Nz")): v.set(d[k])
        self.use_graded.set(d.get("graded",True))
        ms=d.get("mesh") or {}
        self.mesh_type.set(ms.get("type","rect")); self.prism_ax.set(ms.get("prism_axis","auto"))
        self.amr_nodes={k:[float(v) for v in vals] for k,vals in (ms.get("amr") or {}).items() if k in ("x","y","z")}
        self._amr_label()
        self.dense_x,self.dense_y,self.dense_z=d.get("dense_x",[]),d.get("dense_y",[]),d.get("dense_z",[])
        self.regs,self.cons,self.sheets,self.defects=regs,cons,sheets,defects
        sw=d.get("sweep",{})
        self.ivS.set(sw.get("Vs",0.)); self.ivE.set(sw.get("Ve",1.)); self.ivN.set(sw.get("N",11)); self.bv_thr.set(sw.get("bv",0.))
        so=d.get("solver",{})
        for v,k in ((self.T,"T"),(self.max_p,"max_p"),(self.max_c,"max_c"),(self.max_g,"max_g"),
                    (self.tol_p,"tol_p"),(self.tol_c,"tol_c"),(self.m_bgn,"bgn"),(self.m_tau,"tau"),(self.m_vsat,"vsat"),
                    (self.m_qc,"qc"),(self.m_lb,"lombardi"),(self.m_tun,"tunnel"),
                    (self.m_te,"metal_te"),(self.m_sbt,"metal_tunnel"),(self.m_ii,"avalanche")):
            v.set(so[k] if k in so else (False if k in ("qc","lombardi","tunnel","avalanche") else
                                         True if k in ("metal_te","metal_tunnel") else v.get()))
        self.tau0.set(str(so.get("tau0","") or ""))
        self._refresh_groups(tgt=sw.get("tgt"),flt=sw.get("flt") or "(none)",out=sw.get("out") or "auto")
        self._rfr(redraw=False); self._rfc(redraw=False); self._rfd(redraw=False); self._show_sheets()
        self._def_centre(); self._def_sc_whole()
        self.t_name.config(text=self.tname); self.t_cat.config(text=label or f"from {os.path.basename(p)}")
        self.t_doc.configure(state="normal"); self.t_doc.delete("1.0","end")
        base=self.tname.split(" [")[0]
        self.t_doc.insert("end",template_doc(base) if base in TEMPLATES else f"Device loaded from {p}")
        self.t_doc.configure(state="disabled"); self.tmpl_btn.config(text=self.tname[:30])
        self.sol=None; self.ivd=None; self._clear_solution_plots(); self._fill_results(None)
        self._fresh_cuts=True
        view=self.tmeta.get("view")
        if view: self.sec_plane.set(view[0]); self._sec_range(view[1])
        else: self._sec_range(None)
        self._mesh_update()
        if label is None: self.sv.set(f"Opened {p}"); self._log(f"Device opened from {p}")

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

    @staticmethod
    def _core_info():
        st=semibuild.read_stamp()
        if not st or os.environ.get("SEMISIM_SO"): return ""
        return (f" (built with {os.path.basename(st.get('cc','?'))}, "
                f"{'OpenMP' if st.get('openmp') else 'single-threaded: no OpenMP'}"
                f"{', tuned for this CPU' if st.get('native') else ''})")

    def _about(self):
        top=self._text_window(f"About SemiSim 3D {VERSION}",720,560)
        fr,t=self._styled_text(top); fr.pack(fill=tk.BOTH,expand=True,padx=6,pady=6)
        body=f"""
# SemiSim 3D {VERSION}
3-D drift-diffusion semiconductor device simulator: a C core (OpenMP) with a Python/Tk interface.
## Status: test version
- first public version, still being tested: the physics is checked against textbook theory (validate.py), not yet against measurements - check important results independently
- project page: https://github.com/kotmariusz/Semiconductor_devices_simulator
- MIT License: free to use, copy, change and share for any purpose; keep the copyright notice; no warranty
- sources of the models and parameters: Help → Sources (docs/SOURCES.md)
## This session
- core library: `{_SO}`{self._core_info()}
- Python {platform.python_version()}, Tcl/Tk {self.tk.call("info","patchlevel")} ({self.tk.call("tk","windowingsystem")}), NumPy {np.__version__}, Matplotlib {matplotlib.__version__}
- OpenMP threads: {os.environ.get('OMP_NUM_THREADS','?')} (of {self.ncpu} available)
- memory: about {BYTES_PER_NODE} bytes per mesh node; Anderson depth {api_aa() or '-'}
- settings file: `{SETTINGS_FILE}`
## Transport core
- box method on rectangular (tensor) meshes or triangular-prism meshes (Delaunay cross-section x 1-D grid)
- Newton-Poisson (multigrid CG) + Scharfetter-Gummel continuity (multigrid BiCGSTAB)
- Gummel map on (ψ, φn, φp) with Anderson mixing, floating-region balance
- adaptive bias continuation with a sweep predictor, adaptive mesh refinement (rectangular mesh)
- cooperative stop, Kirchhoff-exact terminal currents with resolution floor
## Physics
- heterojunctions, oxide gates (tox, Nox, layer stacks), SiO2/Al2O3/HfO2/Si3N4, buried electrodes, interface sheet charge
- metal contacts: measured/pinned Schottky barriers, thermionic emission with image-force lowering, barrier tunnelling
- Masetti (Si) / Caughey-Thomas (other materials) mobility with material temperature exponents, velocity saturation with Canali's temperature laws (Si), Lombardi surface mobility, SRH τ(N) + Auger + radiative, Klaassen band-gap narrowing, impact ionisation
- defects: deep levels (SRH + trapped charge), interface traps, fixed charge, voids and particles, dislocations, radiation damage
- quantum confinement (MLDA), gate tunnelling (Tsu-Esaki/WKB: direct + Fowler-Nordheim)
## Not modelled
- band-to-band tunnelling, optical generation, subband quantisation, incomplete ionisation, Fermi-Dirac statistics
- negative differential mobility (GaAs, InP, GaN), thermionic heterointerfaces, ballistic transport, self-heating, transient trapping
## Build
- the core is built from `semiconductor_3d.c` on the first start (and after an update of it) by `semibuild.py`: gcc or clang, -O3 -ffast-math, OpenMP (on macOS through Homebrew's libomp); `python3 semibuild.py` rebuilds it by hand
- `semimesh.py` (the triangular-prism mesher), `semidefects.py` (defects) and `semibuild.py` must sit next to main.py
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
