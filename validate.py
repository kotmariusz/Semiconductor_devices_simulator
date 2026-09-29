#!/usr/bin/env python3
"""
SemiSim 3D v7.4 - physics validation battery.

Every section builds a device directly on the C core, solves it and compares
the result with closed-form semiconductor theory.  Each check prints PASS or
FAIL with the tolerance used; the script exits non-zero if anything fails.

    python3 validate.py                 # sections 1-8, 10-13 (about 5-10 minutes)
    python3 validate.py --only 3,10     # selected sections
    python3 validate.py --templates     # + section 9: every GUI template solves

Needs semiconductor_core.so, main.py and semimesh.py in the same folder
(main.py supplies the ctypes bindings and the mesh helpers; no window is
opened).
"""
import os, sys, time, ctypes, argparse
import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import main as M  # noqa: E402

L = M._lib
D3 = ctypes.c_double
_I, _D = ctypes.c_int, ctypes.c_double
for fn, res, args in [
    ("api3_solve", _I, []), ("api3_resolve", _I, []),
    ("api3_set_models", None, [_I]), ("api3_set_contact_V", None, [_I, _D]),
    ("api3_get_contact_I", _D, [_I]), ("api3_get_contact_A", _D, [_I]),
    ("api3_get_contact_Iraw", _D, [_I]), ("api3_get_contact_Ifloor", _D, [_I]),
    ("api3_mat_ni", _D, [_I]), ("api3_get_stat", ctypes.c_long, [_I]),
    ("api3_add_contact_ex", _I, [_I] * 7 + [_I, _D, _D, _D, _D]),
    ("api3_set_gate_stack", _I, [_I, _I, ctypes.POINTER(_I), ctypes.POINTER(_D)]),
    ("api3_get_ntp", _I, [_I]), ("api3_get_Cref", _D, []),
    ("api3_tun_eval", _D, [_I, _I, ctypes.POINTER(_I), ctypes.POINTER(_D), _D, _D, _D, _D, _D, ctypes.POINTER(_D)]),
    ("api3_mlda_factor", _D, [_I, _D, _D, _I]),
    ("api3_lombardi", _D, [_I, _D, _D, _D, _D, _I]),
]:
    f = getattr(L, fn); f.restype = res; f.argtypes = args

Q = 1.602176634e-19; KB = 1.380649e-23; EPS0 = 8.8541878128e-12
VT = KB * 300.0 / Q
SOLV = (2000, 1000, 500, 1e-10, 1e-8, 1.0, 1.0)   # tight tolerances for validation

RESULTS = []


def check(sec, name, value, ref, tol, mode="rel", info=""):
    """mode rel: |value/ref-1|<=tol; abs: |value-ref|<=tol; range: ref<=value<=tol."""
    if mode == "rel":
        err = value / ref - 1.0; ok = abs(err) <= tol
        txt = f"{value:.6g} vs {ref:.6g}  ({err*100:+.3f} %, tol {tol*100:g} %)"
    elif mode == "abs":
        err = value - ref; ok = abs(err) <= tol
        txt = f"{value:.6g} vs {ref:.6g}  (diff {err:+.2e}, tol {tol:g})"
    else:
        ok = ref <= value <= tol
        txt = f"{value:.6g} in [{ref:g}, {tol:g}]"
    RESULTS.append((sec, name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {txt} {info}")
    return ok


# ------------------------------------------------------------------ helpers
def stack1d(layers, xs_m, contacts, ny=3, width=50e-9, solver=SOLV, models=0):
    """1-D bar along x. layers: (x0_nm, x1_nm, material, Nd_cm3, Na_cm3);
    contacts: (face, bc, V, phi_m). Returns (xs, ys, zs) in metres."""
    xs = np.asarray(xs_m, float); Nx = len(xs)
    ys = np.linspace(0, width, ny)
    M.api_init()
    assert M.api_gridc(Nx, ny, ny, (D3 * Nx)(*xs), (D3 * ny)(*ys), (D3 * ny)(*ys))
    M.api_solver(*solver)
    mids = {}
    for (_, _, mat, _, _) in layers:
        if mat not in mids:
            mids[mat] = M.api_addmat(M.MATS[mat], 300.0)
    M.api_region(0, 0, 0, Nx - 1, ny - 1, ny - 1, mids[layers[0][2]], 0.0, 0.0)
    for (x0, x1, mat, Nd, Na) in layers:
        i0, i1 = M.span_idx(xs, x0 * 1e-9, x1 * 1e-9)
        M.api_region(i0, 0, 0, i1, ny - 1, ny - 1, mids[mat], Nd * 1e6, Na * 1e6)
    M.api_clrcon()
    for (face, bc, V, pm) in contacts:
        M.api_addcon(face, 0, ny - 1, 0, ny - 1, bc, V, pm)
    L.api3_set_models(models)
    return xs, ys, ys


def field(which, shape):
    Nz, Ny, Nx = shape; N = Nx * Ny * Nz
    buf = (D3 * N)(); M.api_bulk(which, buf)
    return np.frombuffer(buf, dtype=np.float64, count=N).copy().reshape(Nz, Ny, Nx)


def line(which, xs, ys, zs):
    a = field(which, (len(zs), len(ys), len(xs)))
    return a[len(zs) // 2, len(ys) // 2, :]


PHI, NN_, PP, EC, EV, MUN, MUP = 0, 1, 2, 3, 4, 10, 11
JNX, JPX, NIE = 20, 23, 26


def geo_grid(L_m, h0, growth, hmax):
    x = [0.0]; h = h0
    while x[-1] < L_m:
        x.append(x[-1] + h); h = min(h * growth, hmax)
    x = np.array(x); x *= L_m / x[-1]
    return x


# ------------------------------------------------------------------ sections
def s1_equilibrium():
    print("1. Equilibrium p-n junction (Si): built-in potential and mass action")
    for Na, Nd in [(1e17, 1e17), (1e16, 1e18), (1e15, 1e15)]:
        xs = np.linspace(0, 2e-6, 201)
        stack1d([(0, 1000, "Si", 0, Na), (1000, 2000, "Si", Nd, 0)], xs,
                [(0, 0, 0.0, 0.0), (1, 0, 0.0, 0.0)])
        ok = L.api3_solve()
        phi = line(PHI, xs, xs[:3], xs[:3]); n = line(NN_, xs, xs[:3], xs[:3]); p = line(PP, xs, xs[:3], xs[:3])
        ni = L.api3_mat_ni(0)
        Vbi = VT * np.log(Na * Nd * 1e12 / ni ** 2)
        check(1, f"Vbi Na={Na:.0e} Nd={Nd:.0e}", phi[-1] - phi[0], Vbi, 1e-4, "abs",
              "" if ok else "(not converged)")
        check(1, f"max|np/ni^2-1| Na={Na:.0e} Nd={Nd:.0e}", float(np.max(np.abs(n * p / ni ** 2 - 1))), 0., 1e-9, "range")


def s2_reverse():
    print("2. Reverse-biased abrupt junction: depletion width vs depletion approximation")
    EPS = 11.7 * EPS0; Na = Nd = 1e16; Lnm = 12000.
    xs = np.linspace(0, Lnm * 1e-9, 481)
    for VR in [1, 10, 30, 60]:
        stack1d([(0, Lnm / 2, "Si", 0, Na), (Lnm / 2, Lnm, "Si", Nd, 0)], xs,
                [(0, 0, -VR, 0.0), (1, 0, 0.0, 0.0)])
        ok = L.api3_solve(); ni = L.api3_mat_ni(0)
        n = line(NN_, xs, xs[:3], xs[:3]); p = line(PP, xs, xs[:3], xs[:3])
        half = xs <= xs[-1] / 2
        W = np.trapezoid(np.where(half, 1 - p / (Na * 1e6), 0), xs) + np.trapezoid(np.where(~half, 1 - n / (Nd * 1e6), 0), xs)
        Vbi = VT * np.log(Na * Nd * 1e12 / ni ** 2)
        Wth = np.sqrt(2 * EPS * (Vbi + VR - 2 * VT) / Q * (1 / (Na * 1e6) + 1 / (Nd * 1e6)))
        check(2, f"W at VR={VR} V", W, Wth, 2e-3, "rel", "" if ok else "(not converged)")


def s3_forward():
    print("3. Forward long-base p+/n diode: current vs Shockley (coth) theory")
    EPS = 11.7 * EPS0; Na, Nd = 1e20, 1e16; Wp, Wn = 300., 20000.; Lnm = Wp + Wn
    xs = M.graded_grid(Lnm, 700, [Wp])
    stack1d([(0, Wp, "Si", 0, Na), (Wp, Lnm, "Si", Nd, 0)], xs,
            [(0, 0, 0.0, 0.0), (1, 0, 0.0, 0.0)])
    L.api3_solve(); ni = L.api3_mat_ni(0); Vbi = VT * np.log(Na * Nd * 1e12 / ni ** 2)
    tau = 1e-5
    for V in [0.35, 0.40, 0.45, 0.50, 0.55]:
        L.api3_set_contact_V(0, V); ok = L.api3_resolve()
        mp = line(MUP, xs, xs[:3], xs[:3]); mn = line(MUN, xs, xs[:3], xs[:3])
        mup = np.median(mp[(xs > Wp * 1e-9 + 1e-6) & (xs < xs[-1] - 1e-6)])
        mun = np.median(mn[xs < Wp * 1e-9 - 0.05e-6])
        Dp, Dn = mup * VT, mun * VT; Lp = np.sqrt(Dp * tau)
        xn = np.sqrt(2 * EPS * max(Vbi - V - 2 * VT, 1e-3) / (Q * Nd * 1e6)); W = Wn * 1e-9 - xn
        Jth = (Q * ni ** 2 * Dp / (Nd * 1e6 * Lp) / np.tanh(W / Lp) + Q * ni ** 2 * Dn / (Na * 1e6 * Wp * 1e-9)) * np.expm1(V / VT)
        J = L.api3_get_contact_I(1) / L.api3_get_contact_A(1) * -1.0
        if V < 0.44:
            print(f"  [info] V={V:.2f}: J/J_diffusion = {J/Jth:.4f}  (excess = SRH recombination in the depletion region, expected)")
        else:
            check(3, f"J at V={V:.2f} V", J, Jth, 1e-2, "rel", "" if ok else "(not converged)")


def s4_sic_blocking():
    print("4. 4H-SiC p+/n-/n+ blocking junction to 2 kV: peak field vs depletion approximation")
    EPS = 9.7 * EPS0; Na, Nd = 1e19, 5e15; Wp, Wd, Wn = 500., 25000., 500.; Lnm = Wp + Wd + Wn
    xs = M.graded_grid(Lnm, 420, [Wp, Wp + Wd])
    stack1d([(0, Wp, "4H-SiC", 0, Na), (Wp, Wp + Wd, "4H-SiC", Nd, 0), (Wp + Wd, Lnm, "4H-SiC", 1e19, 0)], xs,
            [(0, 0, 0.0, 0.0), (1, 0, 0.0, 0.0)], solver=(2000, 1000, 500, 1e-8, 1e-6, 1.0, 1.0))
    L.api3_solve(); ni = L.api3_mat_ni(0); Vbi = VT * np.log(Na * Nd * 1e12 / ni ** 2)
    for VR in [100., 1000., 2000.]:
        L.api3_set_contact_V(0, -VR); t = time.time(); ok = L.api3_resolve()
        phi = line(PHI, xs, xs[:3], xs[:3])
        m = (xs > Wp * 1e-9) & (xs < (Wp + Wd) * 1e-9)
        Wth = np.sqrt(2 * EPS * (Vbi + VR - 2 * VT) / (Q * Nd * 1e6))
        Emax = np.max(np.abs(np.gradient(phi, xs)[m])); Eth = Q * Nd * 1e6 * Wth / EPS
        check(4, f"Emax at VR={VR:.0f} V", Emax, Eth, 5e-3, "rel",
              f"[{time.time()-t:.1f} s]" + ("" if ok else " (not converged)"))


def s5_hetero():
    print("5. Abrupt heterojunctions (Anderson rule): built-in potential and band offsets")
    tab = {"GaAs": (1.519, 5.405e-4, 204., 4.07, 4.37e23, 8.68e24), "AlGaAs": (1.882, 5.58e-4, 295., 3.827, 7.2e23, 1.05e25),
           "GaN": (3.507, 9.09e-4, 830., 4.1, 2.3e24, 4.6e25), "AlGaN": (3.990, 9.09e-4, 830., 3.76, 3.2e24, 4.6e25)}

    def mat(name):
        E0, a, b, chi, Nc, Nv = tab[name]; return E0 - a * 300. ** 2 / (300. + b), chi, Nc, Nv
    for (mn, Nd, mp, Na) in [("AlGaAs", 1e17, "GaAs", 1e17), ("AlGaN", 1e18, "GaN", 1e17)]:
        xs = np.linspace(0, 2e-6, 401)
        stack1d([(0, 1000, mn, Nd, 0), (1000, 2000, mp, 0, Na)], xs, [(0, 0, 0.0, 0.0), (1, 0, 0.0, 0.0)])
        L.api3_solve()
        phi = line(PHI, xs, xs[:3], xs[:3]); Ec = line(EC, xs, xs[:3], xs[:3])
        Eg_n, chi_n, Nc_n, _ = mat(mn); Eg_p, chi_p, _, Nv_p = mat(mp)
        Wn = chi_n + VT * np.log(Nc_n / (Nd * 1e6)); Wp = chi_p + Eg_p - VT * np.log(Nv_p / (Na * 1e6))
        i = np.searchsorted(xs, 1000e-9)
        check(5, f"Vbi {mn}/{mp}", phi[0] - phi[-1], Wp - Wn, 1e-4, "abs")
        check(5, f"dEc {mn}/{mp} = dchi", Ec[i] - Ec[i - 1] + (phi[i] - phi[i - 1]), chi_n - chi_p, 1e-3, "abs")


MOS = {}


def s6_mos():
    print("6. MOS capacitor (p-Si 1e17, 10 nm SiO2, n+ poly gate): flat band, threshold, inversion charge")
    EPS = 11.7 * EPS0; EOX = 3.9 * EPS0
    Na = 1e17 * 1e6; tox = 10e-9; phim = 4.05
    x = geo_grid(1.5e-6, 0.1e-9, 1.06, 20e-9); Nx = len(x); ys = np.linspace(0, 20e-9, 3)
    M.api_init(); M.api_gridc(Nx, 3, 3, (D3 * Nx)(*x), (D3 * 3)(*ys), (D3 * 3)(*ys)); M.api_solver(*SOLV)
    L.api3_set_models(0)
    si = M.api_addmat(0, 300.0); M.api_region(0, 0, 0, Nx - 1, 2, 2, si, 0.0, Na)
    M.api_clrcon()
    g = L.api3_add_contact_ex(0, 0, 2, 0, 2, 0, 0, 2, 0.0, phim, tox, 0.0)
    L.api3_add_contact_ex(1, 0, 2, 0, 2, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    L.api3_solve(); ni = L.api3_mat_ni(0)
    Eg = 1.170 - 4.73e-4 * 300 ** 2 / 936.; chi = 4.05; Nv = 3.10e25
    Ws = chi + Eg - VT * np.log(Nv / Na); phiF = VT * np.log(Na / ni); Cox = EOX / tox
    VFB = phim - Ws
    VTth = VFB + 2 * phiF + np.sqrt(2 * EPS * Q * Na * 2 * phiF) / Cox
    res = {}
    for VG in [VFB, VTth, VTth + 1.0, VTth + 1.5]:
        L.api3_set_contact_V(g, VG); ok = L.api3_resolve()
        phi = line(PHI, x, ys, ys); n = line(NN_, x, ys, ys)
        res[VG] = (phi[0] - phi[-1], Q * np.trapezoid(n, x), ok)
    check(6, "surface potential at V_FB", res[VFB][0], 0.0, 1e-3, "abs")
    check(6, "surface potential at V_T(theory) = 2 phi_F", res[VTth][0], 2 * phiF, 5e-3, "abs")
    slope = (res[VTth + 1.5][1] - res[VTth + 1.0][1]) / 0.5
    check(6, "dQinv/dVG / Cox in strong inversion (<1: inversion-layer capacitance)", slope / Cox, 0.90, 0.99, "range")
    MOS.update(VG=VTth + 1.0, Qinv=res[VTth + 1.0][1])


def s7_nmos():
    print("7. Long-channel NMOS (L = 1 um): linear-region drain current vs mu*Qinv*VDS/L")
    if "Qinv" not in MOS:
        s6_mos()
    Lx, Ly, W = 2000e-9, 1000e-9, 100e-9
    xs = M.graded_grid(2000., 90, [500., 1500.])
    ys = geo_grid(Ly, 0.2e-9, 1.12, 50e-9); zs = np.linspace(0, W, 3)
    Nx, Ny, Nz = len(xs), len(ys), len(zs)
    M.api_init(); M.api_gridc(Nx, Ny, Nz, (D3 * Nx)(*xs), (D3 * Ny)(*ys), (D3 * Nz)(*zs)); M.api_solver(*SOLV)
    L.api3_set_models(0)
    si = M.api_addmat(0, 300.0)
    M.api_region(0, 0, 0, Nx - 1, Ny - 1, Nz - 1, si, 0.0, 1e23)
    for (a, b) in [(0, 500e-9), (1500e-9, Lx)]:
        i0, i1 = M.span_idx(xs, a, b); j0, j1 = M.span_idx(ys, 0, 200e-9)
        M.api_region(i0, j0, 0, i1, j1, Nz - 1, si, 1e26, 0.0)
    M.api_clrcon()
    iS = M.span_idx(xs, 0, 450e-9); iD = M.span_idx(xs, 1550e-9, Lx); iG = M.span_idx(xs, 500e-9, 1500e-9)
    cS = L.api3_add_contact_ex(2, iS[0], iS[1], 0, Nz - 1, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    cD = L.api3_add_contact_ex(2, iD[0], iD[1], 0, Nz - 1, 0, 0, 0, 0.05, 0.0, 0.0, 0.0)
    cG = L.api3_add_contact_ex(2, iG[0], iG[1], 0, Nz - 1, 0, 0, 2, MOS["VG"], 4.05, 10e-9, 0.0)
    cB = L.api3_add_contact_ex(3, 0, Nx - 1, 0, Nz - 1, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    ok = L.api3_solve()
    ID = L.api3_get_contact_I(cD)
    mu = float(field(MUN, (Nz, Ny, Nx))[1, 0, Nx // 2])     # low-field channel mobility (Arora, 1e17)
    th = mu * MOS["Qinv"] * 0.05 / 1e-6
    check(7, "ID/W (linear, VDS=50 mV)", ID / W, th, 3e-2, "rel",
          "(theory ignores the channel-charge drop along the channel)" + ("" if ok else " (not converged)"))
    raw = sum(L.api3_get_contact_Iraw(c) for c in (cS, cD, cG, cB))
    check(7, "Kirchhoff |sum I_raw|/ID", abs(raw / ID), 0., 1e-6, "range")


def s8_bjt():
    print("8. 2-D n+pn Si BJT (BGN + tau(N)): electron current in the base vs the Gummel-number integral")
    Lx, Ly, Lz = 2000., 2000., 100.
    xs = M.graded_grid(Lx, 110, [200., 500., 1600.]); ys = M.graded_grid(Ly, 56, [1000., 1400.]); zs = M.uniform_grid(Lz, 3)
    Nx, Ny, Nz = len(xs), len(ys), len(zs)
    M.api_init(); assert M.api_gridc(Nx, Ny, Nz, (D3 * Nx)(*xs), (D3 * Ny)(*ys), (D3 * Nz)(*zs))
    M.api_solver(*SOLV); L.api3_set_models(3)
    si = M.api_addmat(0, 300.0)

    def reg(x0, x1, y0, y1, Nd, Na):
        i0, i1 = M.span_idx(xs, x0 * 1e-9, x1 * 1e-9); j0, j1 = M.span_idx(ys, y0 * 1e-9, y1 * 1e-9)
        M.api_region(i0, j0, 0, i1, j1, Nz - 1, si, Nd * 1e6, Na * 1e6)
    reg(0, Lx, 0, Ly, 1e16, 0); reg(1600, Lx, 0, Ly, 1e19, 0); reg(0, 500, 0, Ly, 0, 1e17)
    reg(0, 200, 0, 1000, 2e19, 0); reg(0, 200, 1400, Ly, 0, 1e19)
    M.api_clrcon()
    jE = M.span_idx(ys, 0, 1000e-9); jB = M.span_idx(ys, 1400e-9, Ly * 1e-9)
    M.api_addcon(0, jE[0], jE[1], 0, Nz - 1, 0, 0.0, 0.0)
    M.api_addcon(0, jB[0], jB[1], 0, Nz - 1, 0, 0.0, 0.0)
    M.api_addcon(1, 0, Ny - 1, 0, Nz - 1, 0, 2.0, 0.0)
    L.api3_solve(); ni = L.api3_mat_ni(0)
    jy = np.searchsorted(ys, 250e-9)
    for VBE in [0.45, 0.55, 0.65, 0.70]:
        L.api3_set_contact_V(1, VBE); ok = L.api3_resolve()
        sh = (Nz, Ny, Nx)
        p = field(PP, sh)[1, jy, :]; nie = field(NIE, sh)[1, jy, :]
        Dn = field(MUN, sh)[1, jy, :] * VT; Jnx = field(JNX, sh)[1, jy, :]
        m = (xs >= 200e-9) & (xs <= 500e-9)
        G = np.trapezoid(p[m] / (Dn[m] * (nie[m] / ni) ** 2), xs[m])
        Jth = Q * ni ** 2 * np.expm1(VBE / VT) / G
        Jsim = -Jnx[np.argmin(np.abs(xs - 350e-9))]
        check(8, f"Jn(base)/Gummel at VBE={VBE:.2f}", Jsim, Jth, 1e-2, "rel", "" if ok else "(not converged)")
        IC = L.api3_get_contact_I(2)
        raw = sum(L.api3_get_contact_Iraw(c) for c in range(3))
        check(8, f"Kirchhoff |sum I_raw|/IC at VBE={VBE:.2f}", abs(raw / IC), 0., 1e-4 if VBE < 0.5 else 1e-6, "range")
        if VBE == 0.65:
            IB = L.api3_get_contact_I(1)
            print(f"  [info] beta = IC/IB = {IC/IB:.1f} (with band-gap narrowing in the 2e19 emitter)")


def s9_templates():
    print("9. GUI templates: every built-in device solves at its default bias")
    for name, fn in M.TEMPLATES.items():
        t = fn(); t0 = time.time()
        try:
            M.setup_device(t, models=t.get("models", 7))
            ok = M.api_solve()
        except Exception as e:  # noqa: BLE001
            check(9, name, 0, 1, 0, "range", f"exception {e}"); continue
        dt = time.time() - t0
        car = list(range(len(t["cons"])))            # gates carry the tunnelling current (0 without it)
        raw = sum(L.api3_get_contact_Iraw(k) for k in car)
        Imax = max([abs(L.api3_get_contact_Iraw(k)) for k in car] + [1e-300])
        fl = max([L.api3_get_contact_Ifloor(k) for k in car] + [0.])
        kir = abs(raw) / Imax if Imax > 10 * fl else 0.0
        check(9, f"{name} converged", float(ok), 1.0, 1.0, "range",
              f"[{dt:.1f} s, {L.api3_get_stat(3)} Gummel its, Kirchhoff {kir:.1e}]")


def s10_cmos_finfet():
    """CMOS inverter (floating-output transfer curve) and bulk FinFET."""
    print("10. CMOS inverter and FinFET: thresholds, transfer curve, subthreshold swing, DIBL")
    t = M.TEMPLATES["CMOS Inverter (0.25um)"](); cons = t["cons"]
    M.setup_device(t, models=7)          # classical physics: compared with long-channel theory
    vin, out, vdd = M.group_of(cons, "Vin"), M.group_of(cons, "Out"), M.group_of(cons, "VDD")
    # long-channel theory: V_th = V_FB + 2 phi_F + sqrt(2 q eps Na 2phi_F)/Cox, Na = 5e17, tox 5 nm
    ni = L.api3_mat_ni(0) * 1e-6; phiF = VT * np.log(5e17 / ni)
    Cox = 3.9 * EPS0 / 5e-9; Qd = np.sqrt(2 * Q * 11.7 * EPS0 * 5e23 * 2 * phiF)
    Eg = 1.12
    VtN = 4.05 - (4.05 + Eg / 2 + phiF) + 2 * phiF + Qd / Cox
    VtP = 5.17 - (4.05 + Eg / 2 - phiF) - 2 * phiF - Qd / Cox
    va = [c.V for c in cons]
    for c in out: va[c] = 0.1
    R = M.run_sweep(cons, vin, 0., 1.2, 13, vapp=va)
    I = R["I"][:, out[0]]; V = R["V"]; gm = np.gradient(I, V); j = int(np.argmax(gm))
    check(10, "NMOS V_th (max-g_m, V_DS = 0.1 V) vs long-channel theory", V[j] - I[j] / gm[j] - 0.05, VtN, 0.06, "abs")
    for c in out: va[c] = 2.4
    R = M.run_sweep(cons, vin, 2.5, 1.3, 13, vapp=va)
    I = -R["I"][:, out[1]]; V = R["V"]; gm = -np.gradient(I, V); j = int(np.argmax(gm))
    check(10, "PMOS V_th (max-g_m, V_SD = 0.1 V) vs long-channel theory", V[j] + I[j] / gm[j] - 2.5 + 0.05, VtP, 0.06, "abs")
    va = [c.V for c in cons]
    R = M.run_sweep(cons, vin, 0., 2.5, 11, flt=out, vapp=va, refine=False)
    Vi, Vo = R["V"], R["Vf"]
    Iout = np.abs(R["I"][:, out].sum(axis=1)); Isup = np.abs(R["I"][:, vdd[0]])
    Fl = R["F"][:, out].sum(axis=1)
    check(10, "VTC: every point solved (floating output, I = 0)", float(R["conv"].all()), 1., 1., "range")
    check(10, "VTC: V_OH at V_in = 0 (V)", Vo[0], 2.499, 2.5 + 1e-9, "range")
    check(10, "VTC: V_OL at V_in = VDD (V)", Vo[-1], -1e-9, 1e-3, "range")
    bal = np.max(Iout / np.maximum(1e-3 * Isup + 3 * Fl, 1e-300))
    check(10, "VTC: output current balance |I_out| / tolerance", bal, 0., 1.0 + 1e-9, "range")
    d_ = Vo - Vi; k = np.where(np.diff(np.sign(d_)) != 0)[0]
    VM = Vi[k[0]] + (Vi[k[0] + 1] - Vi[k[0]]) * d_[k[0]] / (d_[k[0]] - d_[k[0] + 1]) if len(k) else np.nan
    # square law with the channel mobilities of 5e17 Si (mu_n/mu_p ~ 2.5): V_M ~ 1.07 V
    check(10, "VTC: switching threshold V_M (V), square-law estimate 1.07 V", VM, 0.95, 1.25, "range")
    ok_mono = bool(np.all(np.diff(Vo) <= 1e-6))
    check(10, "VTC: V_out non-increasing in V_in", float(ok_mono), 1., 1., "range")

    t = M.TEMPLATES["FinFET (bulk tri-gate)"](); cons = t["cons"]
    M.setup_device(t, models=t.get("models", 7))     # as the GUI runs it (MLDA, Lombardi, tunnelling)
    gate, dr = M.group_of(cons, "Gate"), M.group_of(cons, "Drain")[0]
    va = [c.V for c in cons]
    Rs = M.run_sweep(cons, gate, 0., 0.8, 17, vapp=va)
    va[dr] = 0.05
    Rl = M.run_sweep(cons, gate, 0.8, 0., 17, vapp=va)
    Is, Vs_ = np.abs(Rs["I"][:, dr]), Rs["V"]
    Il, Vl = np.abs(Rl["I"][:, dr])[::-1], Rl["V"][::-1]
    lg = np.log10(Is); ss = np.diff(Vs_) / np.diff(lg) * 1e3
    SS = float(np.min(ss[Is[1:] < 1e-3 * Is.max()]))
    check(10, "FinFET SS at V_D = 0.8 V (mV/dec): >= thermionic limit 59.6", SS, 59.5, 75., "range")
    Vat = lambda V_, I_, Ic: float(np.interp(np.log10(Ic), np.log10(I_), V_))
    DIBL = (Vat(Vl, Il, 1e-9) - Vat(Vs_, Is, 1e-9)) / 0.75 * 1e3
    check(10, "FinFET DIBL (mV/V) at I_D = 1 nA", DIBL, 0., 60., "range")
    check(10, "FinFET I_on/I_off (V_D = 0.8 V)", Is[-1] / Is[0], 1e5, 1e9, "range")
    Ig = abs(Rs["I"][-1, gate].sum())
    check(10, "FinFET HfO2-stack gate leakage I_G/I_D at V_G = V_D = 0.8 V", Ig / Is[-1], 1e-10, 1e-4, "range")
    print(f"  [info] FinFET: I_on = {Is[-1]:.3e} A/fin, I_off = {Is[0]:.3e} A/fin, SS = {SS:.1f} mV/dec, "
          f"DIBL = {DIBL:.1f} mV/V, I_G(0.8 V) = {Ig:.2e} A")


# ------------------------------------------------------------ nanoscale MOS
def nmos2d(stack, VG, models, Na=1e18, pm=4.15, Lg=100e-9, tail=None):
    """Planar n-MOSFET cross-section with an oxide (Robin) gate of the given
    stack [(material, m)], grounded n+ source/drain. Returns (ok, J_G A/cm^2,
    contact ids, xs, ys)."""
    Lx = Lg + 200e-9; Ly = 150e-9
    xs = M.graded_grid(Lx * 1e9, 44, [100., 100. + Lg * 1e9])
    ys = geo_grid(Ly, 0.1e-9, 1.12, 10e-9); zs = np.linspace(0, 20e-9, 3)
    Nx, Ny, Nz = len(xs), len(ys), len(zs)
    M.api_init(); M.api_gridc(Nx, Ny, Nz, (D3 * Nx)(*xs), (D3 * Ny)(*ys), (D3 * Nz)(*zs)); M.api_solver(*SOLV)
    L.api3_set_models(models)
    si = M.api_addmat(0, 300.0)
    M.api_region(0, 0, 0, Nx - 1, Ny - 1, Nz - 1, si, 0.0, Na * 1e6)
    for (a, b) in [(0, 100e-9), (100e-9 + Lg, Lx)]:
        i0, i1 = M.span_idx(xs, a, b); j0, j1 = M.span_idx(ys, 0, 50e-9)
        M.api_region(i0, j0, 0, i1, j1, Nz - 1, si, 1e26, 0.0)
    ids = [M.api_addmat(M.MATS[m], 300.0) for m, _ in stack]
    eot = sum(t * 3.9 / M.EPS_R[m] for m, t in stack)
    M.api_clrcon()
    iS = M.span_idx(xs, 0, 90e-9); iD = M.span_idx(xs, 110e-9 + Lg, Lx); iG = M.span_idx(xs, 100e-9, 100e-9 + Lg)
    cS = L.api3_add_contact_ex(2, iS[0], iS[1], 0, Nz - 1, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    cD = L.api3_add_contact_ex(2, iD[0], iD[1], 0, Nz - 1, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    cG = L.api3_add_contact_ex(2, iG[0], iG[1], 0, Nz - 1, 0, 0, 2, 0.0, pm, eot, 0.0)
    cB = L.api3_add_contact_ex(3, 0, Nx - 1, 0, Nz - 1, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    n = len(stack)
    L.api3_set_gate_stack(cG, n, (_I * n)(*ids), (D3 * n)(*[t for _, t in stack]))
    L.api3_solve()
    out = []
    for v in (VG if hasattr(VG, "__len__") else [VG]):
        L.api3_set_contact_V(cG, v); ok = L.api3_resolve()
        out.append((ok, L.api3_get_contact_I(cG) / L.api3_get_contact_A(cG) * 1e-4))
        if tail: tail(v, xs, ys)
    return out, (cS, cD, cG, cB), xs, ys


def s11_tunnelling():
    print("11. Gate tunnelling: WKB/Tsu-Esaki kernel, direct tunnelling vs thickness, high-k, Fowler-Nordheim")
    # (a) kernel vs an independent integration: flat 1.2 nm SiO2 barrier (V_ox = 0)
    M.api_init(); si = M.api_addmat(0, 300.0); ox = M.api_addmat(8, 300.0); out = (D3 * 4)()
    L.api3_tun_eval(si, 1, (_I * 1)(ox), (D3 * 1)(1e-9), 0., 0., 0., 0., 4.5, out)   # prepares the gauge
    Cref = L.api3_get_Cref(); chi_s, chi_o, m_ox, m_s = 4.05, 0.95, 0.42, 0.50
    HB, M0_, H = 1.054571817e-34, 9.1093837015e-31, 6.62607015e-34
    def ref(EF, Ec, t_ox, u0, u1):
        """electron flux for a barrier (above E = Ec) linear from u0 at the interface to u1 at
        the metal: WKB exponent by numerical quadrature in x, trapezoid in energy"""
        E = np.linspace(Ec, max(EF, Ec) + 25 * VT, 4001)
        x = np.linspace(0., 1., 4001)
        u = np.maximum((u0 + (u1 - u0) * x)[None, :] - (E - Ec)[:, None], 0.0)
        T = np.exp(-2 * np.trapezoid(np.sqrt(2 * m_ox * M0_ * Q * u), x, axis=1) * t_ox / HB)
        xx = (EF - E) / VT
        sp = np.where(xx > 35, xx + np.log1p(np.exp(-np.abs(xx))), np.log1p(np.exp(np.minimum(xx, 35))))
        return 4 * np.pi * m_s * M0_ * Q * Q * VT / H ** 3 * np.trapezoid(T * sp, E)
    # direct: flat 1.2 nm barrier (V_ox = 0); psi_m = Vg + Cref - pm with pm = Cref
    psi = 0.35; t_ox = 1.2e-9; Vg = psi
    L.api3_tun_eval(si, 1, (_I * 1)(ox), (D3 * 1)(t_ox), psi, 0., 0., Vg, Cref, out)
    Ec = -psi - chi_s; phB = chi_s - chi_o
    check(11, "Tsu-Esaki/WKB kernel, direct (flat barrier): flux s->g vs integration", out[0],
          ref(-Cref, Ec, t_ox, phB, phB), 2e-3, "rel")
    check(11, "Tsu-Esaki/WKB kernel, direct (flat barrier): flux g->s vs integration", out[1],
          ref(-Vg - Cref, Ec, t_ox, phB, phB), 2e-3, "rel")
    # Fowler-Nordheim: 6 nm, 6 V across the oxide - the barrier ends inside the oxide
    t_ox = 6e-9; Vox = 6.0; Vg = psi + Vox
    L.api3_tun_eval(si, 1, (_I * 1)(ox), (D3 * 1)(t_ox), psi, 0., 0., Vg, Cref, out)
    check(11, "Tsu-Esaki/WKB kernel, Fowler-Nordheim (triangular barrier): flux s->g vs integration", out[0],
          ref(-Cref, Ec, t_ox, phB, phB - Vox), 2e-3, "rel")
    # (b) self-consistent n-MOSFET, 1.2 nm SiO2, metal gate 4.15 eV, p 1e18
    res = {}
    for tox in (1.0e-9, 1.2e-9, 1.5e-9):
        r, cid, *_ = nmos2d([("SiO2", tox)], [1.0], 32)
        res[tox] = r[0][1]
        if tox == 1.2e-9:
            raw = sum(L.api3_get_contact_Iraw(c) for c in cid); IG = L.api3_get_contact_I(cid[2])
            flr = sum(L.api3_get_contact_Ifloor(c) for c in cid)
            check(11, "Kirchhoff with gate current: |sum I_raw| / (1e-6 I_G + round-off floors)",
                  abs(raw) / (1e-6 * abs(IG) + flr), 0., 1., "range")
    check(11, "J_G(1 V), 1.2 nm SiO2 (A/cm^2): textbook 1e2-1e3", res[1.2e-9], 30., 3000., "range")
    dec = np.log10(res[1.0e-9] / res[1.5e-9]) / 0.5
    check(11, "direct tunnelling: decades per nm of SiO2 (~5: 10x per 0.2 nm)", dec, 3.5, 6.5, "range")
    r, *_ = nmos2d([("SiO2", 1.2e-9)], [1.0], 40)
    check(11, "MLDA on/off: J_G(1 V) ratio (supply from the wall pressure)", r[0][1] / res[1.2e-9], 0.5, 2.0, "range")
    # (c) high-k stack of the same EOT (1.2 nm): 0.5 nm SiO2 + 3.95 nm HfO2
    r, *_ = nmos2d([("SiO2", 0.5e-9), ("HfO2", 3.95e-9)], [1.0], 32)
    check(11, "SiO2 / (SiO2+HfO2) leakage ratio at equal EOT 1.2 nm, V_G = 1 V", res[1.2e-9] / r[0][1], 1e2, 1e5, "range")
    # (d) Fowler-Nordheim: 8 nm SiO2; oxide field from the solved surface potential
    E, J = [], []
    def tail(v, xs, ys):
        phi = field(PHI, (3, len(ys), len(xs)))
        i = int(np.argmin(np.abs(xs - (xs[0] + xs[-1]) / 2)))
        E.append((v + L.api3_get_Cref() - 4.15 - phi[1, 0, i]) / 8e-9)
    r, *_ = nmos2d([("SiO2", 8e-9)], [6., 7., 8., 9., 10.], 32, Na=5e17, tail=tail)
    J = np.array([x[1] for x in r]); E = np.array(E)
    B = -np.polyfit(1 / E, np.log(J / E ** 2), 1)[0]
    Bth = 4. / 3. * np.sqrt(2 * 0.42 * M0_) * (3.10 * Q) ** 1.5 / (HB * Q)
    # measured Si/SiO2 Fowler-Nordheim plots give B = 2.3-2.5e10 V/m (Lenzlinger & Snow 1969: 2.5e10);
    # the degenerate inversion-layer supply lowers the fitted slope a few % below the barrier value
    check(11, "Fowler-Nordheim plot slope B (V/m), measured Si/SiO2 2.3-2.5e10", B, 2.2e10, 2.6e10, "range",
          f"[barrier value {Bth:.3e}; J = {J[0]:.1e} .. {J[-1]:.1e} A/cm^2 at {E[0]*1e-8:.1f} .. {E[-1]*1e-8:.1f} MV/cm]")


def mosc(tox, Na, VGs, models, pm=4.05):
    """1-D MOS capacitor (Robin gate, SiO2 tox), p-Si Na: per V_G the inversion
    sheet density (cm^-2), its centroid (nm) and the surface potential."""
    x = geo_grid(1.0e-6, 0.05e-9, 1.05, 10e-9); Nx = len(x); ys = np.linspace(0, 20e-9, 3)
    M.api_init(); M.api_gridc(Nx, 3, 3, (D3 * Nx)(*x), (D3 * 3)(*ys), (D3 * 3)(*ys)); M.api_solver(*SOLV)
    L.api3_set_models(models)
    si = M.api_addmat(0, 300.0); M.api_region(0, 0, 0, Nx - 1, 2, 2, si, 0.0, Na * 1e6)
    M.api_clrcon()
    g = L.api3_add_contact_ex(0, 0, 2, 0, 2, 0, 0, 2, 0.0, pm, tox, 0.0)
    L.api3_add_contact_ex(1, 0, 2, 0, 2, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    L.api3_solve(); out = []
    for v in VGs:
        L.api3_set_contact_V(g, v); L.api3_resolve()
        n = line(NN_, x, ys, ys); m = x < 50e-9
        Ns = np.trapezoid(n[m], x[m]); zc = np.trapezoid(n[m] * x[m], x[m]) / Ns
        out.append((Ns * 1e-4, zc * 1e9))
    return np.array(out)


def s12_mlda_lombardi():
    print("12. Quantum confinement (MLDA) and Lombardi surface mobility")
    # (a) the cell-averaged wall factor vs direct integration
    M.api_init(); si = M.api_addmat(0, 300.0)
    lam = lambda m: 1.054571817e-34 / np.sqrt(2 * m * 9.1093837015e-31 * 1.380649e-23 * 300.)
    worst = 0.
    for a, b in [(0., 0.05e-9), (0., 0.5e-9), (0.3e-9, 1.1e-9), (2e-9, 3e-9)]:
        z = np.linspace(a, b, 200001)
        f = lambda m: np.trapezoid(1 - np.exp(-(z / lam(m)) ** 2), z) / (b - a)
        ref = f(0.916) / 3. + 2. * f(0.19) / 3.           # Delta-2 (1/3) + Delta-4 (2/3) valleys
        worst = max(worst, abs(L.api3_mlda_factor(si, a, b, 0) / ref - 1))
    check(12, f"MLDA cell average (Delta-2 {lam(0.916)*1e9:.2f} nm, Delta-4 {lam(0.19)*1e9:.2f} nm) vs integration, "
              f"worst rel. error", worst, 0., 1e-5, "range")
    # (b) MOS capacitor, 2 nm SiO2, p 5e17: threshold shift, dark space, capacitance thickness
    VG = np.linspace(-0.2, 2.2, 49)
    cl = mosc(2e-9, 5e17, VG, 0); qm = mosc(2e-9, 5e17, VG, 8)
    Vat = lambda R, Ns: float(np.interp(np.log10(Ns), np.log10(R[:, 0]), VG))
    dVt = Vat(qm, 1e11) - Vat(cl, 1e11)
    check(12, "MLDA threshold shift at N_inv = 1e11 cm^-2, p 5e17 (V): QM 0.04-0.10", dVt, 0.03, 0.12, "range")
    zq = float(np.interp(np.log10(1e13), np.log10(qm[:, 0]), qm[:, 1]))
    zc = float(np.interp(np.log10(1e13), np.log10(cl[:, 0]), cl[:, 1]))
    # Fang-Howard (Delta-2, T = 0): 3/b = 0.55 nm here; Delta-4 and excited subbands at 300 K push it out
    check(12, "inversion-charge centroid at 1e13 cm^-2 (nm): Schroedinger-Poisson 0.6-1.2", zq, 0.6, 1.2, "range",
          f"[classical {zc:.2f} nm]")
    k = (cl[:, 0] > 5e12) & (cl[:, 0] < 1.5e13); kq = (qm[:, 0] > 5e12) & (qm[:, 0] < 1.5e13)
    Ccl = np.polyfit(VG[k], cl[k, 0], 1)[0] * Q * 1e4; Cqm = np.polyfit(VG[kq], qm[kq, 0], 1)[0] * Q * 1e4
    dCET = 3.9 * EPS0 * (1 / Cqm - 1 / Ccl) * 1e9
    # full Schroedinger-Poisson adds ~0.1-0.3 nm at these densities; MLDA (a near-wall correction)
    # captures the lower end - reported, checked only for the right sign and order
    check(12, "capacitance-equivalent thickness added by quantisation (nm)", dCET, 0.03, 0.4, "range")
    # (c) Lombardi vs the universal mobility curve: long-channel NMOS, 5 nm SiO2, p 1e17
    Lch, W = 1e-6, 100e-9
    xs = M.graded_grid(2000., 70, [500., 1500.]); ys = geo_grid(1e-6, 0.1e-9, 1.12, 50e-9); zs = np.linspace(0, W, 3)
    Nx, Ny, Nz = len(xs), len(ys), len(zs)
    M.api_init(); M.api_gridc(Nx, Ny, Nz, (D3 * Nx)(*xs), (D3 * Ny)(*ys), (D3 * Nz)(*zs)); M.api_solver(*SOLV)
    L.api3_set_models(16)
    si = M.api_addmat(0, 300.0); M.api_region(0, 0, 0, Nx - 1, Ny - 1, Nz - 1, si, 0.0, 1e23)
    for (a, b) in [(0, 500e-9), (1500e-9, 2000e-9)]:
        i0, i1 = M.span_idx(xs, a, b); j0, j1 = M.span_idx(ys, 0, 150e-9)
        M.api_region(i0, j0, 0, i1, j1, Nz - 1, si, 1e26, 0.0)
    M.api_clrcon()
    iS = M.span_idx(xs, 0, 450e-9); iD = M.span_idx(xs, 1550e-9, 2000e-9); iG = M.span_idx(xs, 500e-9, 1500e-9)
    L.api3_add_contact_ex(2, iS[0], iS[1], 0, Nz - 1, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    cD = L.api3_add_contact_ex(2, iD[0], iD[1], 0, Nz - 1, 0, 0, 0, 0.02, 0.0, 0.0, 0.0)
    cG = L.api3_add_contact_ex(2, iG[0], iG[1], 0, Nz - 1, 0, 0, 2, 1.0, 4.05, 5e-9, 0.0)
    L.api3_add_contact_ex(3, 0, Nx - 1, 0, Nz - 1, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)
    L.api3_solve(); Cref = L.api3_get_Cref(); EPS = 11.7 * EPS0; Cox = 3.9 * EPS0 / 5e-9
    i = int(np.argmin(np.abs(xs - 1000e-9))); worst = 0.; info = []
    for VG in (1.2, 1.6, 2.2, 2.8):
        L.api3_set_contact_V(cG, VG); ok = L.api3_resolve()
        sh = (Nz, Ny, Nx); n = field(NN_, sh)[1, :, i]; phi = field(PHI, sh)[1, :, i]
        m = ys < 30e-9; Qi = Q * np.trapezoid(n[m], ys[m])
        Fs = Cox * (VG + Cref - 4.05 - phi[0]) / EPS            # surface field (Gauss at the oxide)
        Eeff = (Fs - 0.5 * Qi / EPS) * 1e-8                     # MV/cm: (Q_dep + Q_inv/2)/eps
        mu = L.api3_get_contact_I(cD) * Lch / (W * Qi * 0.02) * 1e4
        mu_u = 540. / (1 + (Eeff / 0.9) ** 1.85)                 # universal curve (electrons)
        worst = max(worst, abs(mu / mu_u - 1)); info.append(f"{Eeff:.2f} MV/cm: {mu:.0f} vs {mu_u:.0f}")
    check(12, "Lombardi mu_eff vs universal curve 540/(1+(E/0.9 MV/cm)^1.85), worst rel. dev.", worst, 0., 0.25, "range",
          "[" + "; ".join(info) + "]")


# ------------------------------------------------------------------ 13
def _currents(t, models, mesh, thin=True, volt=None):
    """Solve template dict t at its default bias (volt: {label: V}) on the
    given mesh; returns (terminal currents, node count, seconds)."""
    t = dict(t); t["mesh"] = mesh; t["tri_thin"] = thin
    t0 = time.time()
    M.setup_device(t, solver=SOLV, models=models, volt=volt)
    cons = t["cons"]; vapp = [(volt or {}).get(c.label, c.V) for c in cons]
    ok, _ = M.run_point(cons, [], None, vapp)
    I = [L.api3_get_contact_I(k) for k in range(len(cons))]
    return ok, I, M.api_Nx() * M.api_Ny() * M.api_Nz(), time.time() - t0


def _radial_pb(dV, ni, Na, R, tox, rs, planar=False):
    """Band bending (V) at radii rs (nm) of a p-type Si cylinder (radius R)
    in an SiO2 shell (tox) inside a metal, dV = psi_metal - psi_centre,
    Boltzmann statistics, neutral centre: fine-grid Newton on the radial
    Poisson-Boltzmann equation, the oxide as its exact log-potential."""
    eps_s, eps_o = 11.7 * EPS0, 3.9 * EPS0
    p0 = Na / 2 + np.sqrt((Na / 2) ** 2 + ni ** 2); n0 = ni ** 2 / p0
    Rm = R * 1e-9; n = 3001; r = np.linspace(0, Rm, n); h = r[1] - r[0]
    rf = np.r_[0., 0.5 * (r[1:] + r[:-1]), Rm]
    vol = 0.5 * (rf[1:] ** 2 - rf[:-1] ** 2)
    g = eps_s * rf[1:-1] / h
    Cox = eps_o / (Rm * np.log((R + tox) / R))
    u = np.zeros(n); uB = dV / VT; i = np.arange(n - 1)
    for _ in range(300):
        F = np.zeros(n); d0 = np.zeros(n); dp = np.zeros(n - 1); dm = np.zeros(n - 1)
        F[:-1] += g * (u[1:] - u[:-1]); F[1:] += g * (u[:-1] - u[1:])
        d0[:-1] -= g; d0[1:] -= g; dp += g; dm += g
        F[-1] += Cox * Rm * (uB - u[-1]); d0[-1] -= Cox * Rm
        F += Q / VT * (p0 * np.exp(-u) - n0 * np.exp(u) - Na) * vol
        d0 += Q / VT * (-p0 * np.exp(-u) - n0 * np.exp(u)) * vol
        # tridiagonal solve (Thomas)
        b = d0.copy(); rhs = -F.copy()
        for k in range(1, n):
            w = dm[k - 1] / b[k - 1]; b[k] -= w * dp[k - 1]; rhs[k] -= w * rhs[k - 1]
        du = np.zeros(n); du[-1] = rhs[-1] / b[-1]
        for k in range(n - 2, -1, -1): du[k] = (rhs[k] - dp[k] * du[k + 1]) / b[k]
        du = np.clip(du, -2., 2.); u += du
        if np.max(np.abs(du)) < 1e-11: break
    return np.interp(np.asarray(rs) * 1e-9, r, u) * VT, u[0] * VT


def _cylinder_moscap(mesh, Ny, R=30., tox=2., Na=1e19, VGs=(-1.5, -0.5, 0.5, 1.5, 2.5)):
    """Si cylinder + SiO2 shell + wrapped metal gate, body contact at the axis.
    Returns [(VG, surface band bending sim, ref)] and the node count."""
    c0 = R + tox + 3.; Lw = 2 * c0; Lt = (6., Lw, Lw)
    wire = M.SM.Shape.circle("x", c0, c0, R); ox = wire.offset(tox)
    t = dict(Lx=6., Ly=Lw, Lz=Lw, Nx=3, Ny=Ny, Nz=Ny, dense_x=[], dense_y=[], dense_z=[],
             regs=[M.Reg("SiO2", 0, 0, 0, 6, Lw, Lw, "i", 0., "oxide"),
                   M.Reg("Si", 0, 0, 0, 6, Lw, Lw, "p", Na, "p-Si wire", shape=wire)],
             cons=[M.con_box(0, 0, 0, 6, Lw, Lw, Lt, 2, 0., "Gate", 4.05, shape=ox.outside()),
                   M.con_box(0, c0 - 1, c0 - 1, 6, c0 + 1, c0 + 1, Lt, 0, 0., "Body")],
             mesh=mesh, extrude="x")
    M.setup_device(t, solver=SOLV, models=0)
    ni = L.api3_mat_ni(1)                       # materials: 0 SiO2, 1 Si
    out = []
    for k, VG in enumerate(VGs):
        L.api3_set_contact_V(0, VG)
        (L.api3_solve if k == 0 else L.api3_resolve)()
        s = M.pull()
        if mesh == "tri":
            X = s["prism"].node_xyz(); phi = s["raw"]["phi"]; kind = s["raw"]["kind"]
        else:
            Z3, Y3, X3 = np.meshgrid(s["z"], s["y"], s["x"], indexing="ij")
            X = np.c_[X3.ravel(), Y3.ravel(), Z3.ravel()]; phi = s["phi"].ravel(); kind = s["kind"].ravel()
        mid = np.abs(X[:, 0] - 3.) < 1e-6; r = np.hypot(X[:, 1] - c0, X[:, 2] - c0)
        psi_m = np.median(phi[mid & (kind == 2) & (r > R + tox - 0.01) & (r < R + tox + 2)])
        psi_0 = np.median(phi[mid & (r < 0.8)])
        on = mid & (kind == 0) & (np.abs(X[:, 2] - c0) < 1e-6) & (X[:, 1] > c0)
        j = np.nonzero(on)[0][np.argmax(r[on])]            # outermost Si node on the +y radius
        ref, uc = _radial_pb(psi_m - psi_0, ni, Na * 1e6, R, tox, [r[j]])
        out.append((VG, phi[j] - psi_0, ref[0], uc))
    return out, M.api_Nx() * M.api_Ny() * M.api_Nz()


def s13_prism_mesh():
    """Triangular-prism meshes: the general (graph) box-method path, the
    thinned prism mesh against the tensor mesh, round geometry against the
    radial Poisson-Boltzmann solution, and the GAA nanowire template."""
    print("13. Triangular-prism mesh: graph path = tensor path, thinned mesh, round MOS capacitor, GAA FET")
    # 13a: every tensor point kept -> the prism mesh IS the tensor mesh, triangulated
    for name, models, volt in (("Gate leakage: 1.2nm SiO2 NMOS", 47, {"Drain": 0.5}),
                               ("FinFET (bulk tri-gate)", 7, None)):
        t = M.TEMPLATES[name]()
        _, Ir, Nr, tr = _currents(t, models, "rect", volt=volt)
        _, It, Nt, tt = _currents(t, models, "tri", thin=False, volt=volt)
        k = int(np.argmax(np.abs(Ir)))
        check(13, f"{name}: prism mesh with every tensor node, {t['cons'][k].label} current", It[k], Ir[k], 1e-5,
              info=f"[{Nt} nodes; {tr:.1f} s / {tt:.1f} s]")
    # 13b: the thinned prism mesh against the tensor mesh (same physics)
    for name, tol, volt in (("NMOS", 0.01, None), ("Gate leakage: 1.2nm SiO2 NMOS", 0.01, None),
                            ("GaN HEMT", 0.005, None), ("FinFET (bulk tri-gate)", 0.015, None)):
        t = M.TEMPLATES[name](); mdl = t.get("models", 7)
        _, Ir, Nr, tr = _currents(t, mdl, "rect", volt=volt)
        _, It, Nt, tt = _currents(t, mdl, "tri", volt=volt)
        for k in sorted(range(len(Ir)), key=lambda q: -abs(Ir[q]))[:2]:
            if abs(Ir[k]) < 1e-12: continue
            check(13, f"{name}: thinned prism mesh, {t['cons'][k].label} current", It[k], Ir[k], tol,
                  info=f"[{Nt} vs {Nr} nodes; {tt:.1f} s vs {tr:.1f} s]")
    # 13c: round MOS capacitor, prism mesh vs radial Poisson-Boltzmann; staircase for comparison
    tri, Nt = _cylinder_moscap("tri", 49)
    stc, Ns = _cylinder_moscap("rect", 49)
    for (VG, a, b, uc), (_, a2, b2, _) in zip(tri, stc):
        check(13, f"cylinder MOS (R 30 nm, 2 nm SiO2, 1e19) V_G = {VG:+.1f} V: surface band bending",
              a, b, 5e-3, "abs", f"[round prism mesh, {Nt} nodes; staircase {Ns} nodes: {1e3*(a2-b2):+.1f} mV]")
    etri = max(abs(a - b) for _, a, b, _ in tri); estc = max(abs(a - b) for _, a, b, _ in stc)
    check(13, "round prism mesh at least 5x more accurate than the staircase tensor mesh",
          estc / max(etri, 1e-6), 5., 1e9, "range", f"[max error {1e3*etri:.1f} vs {1e3*estc:.1f} mV]")
    check(13, "cylinder: neutral centre in the reference (depletion < R)", max(abs(uc) for *_, uc in tri), 0., 1e-3, "abs")
    # 13d: GAA nanowire FET template: subthreshold swing of a gate-all-around wire
    t = M.TEMPLATES["GAA nanowire FET (round)"](); cons = t["cons"]
    M.setup_device(t, solver=SOLV, models=t["models"])
    vapp = [c.V for c in cons]; I = {}
    for VG in (0.05, 0.15, 0.7):
        vapp[2] = VG; ok, _ = M.run_point(cons, [], None, vapp); I[VG] = L.api3_get_contact_I(1)
    SS = 100. / np.log10(I[0.15] / I[0.05])
    check(13, "GAA nanowire (L_g 20 nm, D 10 nm): subthreshold swing (mV/dec), ~60-65 for a gate all round",
          SS, 59., 70., "range", f"[I_D {I[0.05]:.2e} / {I[0.15]:.2e} A]")
    IG = L.api3_get_contact_I(2)
    check(13, "GAA nanowire: I_on/I_off (0.7 V / 0.05 V)", I[0.7] / I[0.05], 1e3, 1e9, "range",
          f"[I_on {I[0.7]*1e6:.1f} uA, gate leakage {IG:.2e} A]")


SECTIONS = {1: s1_equilibrium, 2: s2_reverse, 3: s3_forward, 4: s4_sic_blocking, 5: s5_hetero,
            6: s6_mos, 7: s7_nmos, 8: s8_bjt, 9: s9_templates, 10: s10_cmos_finfet,
            11: s11_tunnelling, 12: s12_mlda_lombardi, 13: s13_prism_mesh}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="", help="comma-separated section numbers")
    ap.add_argument("--templates", action="store_true", help="also run section 9 (all GUI templates)")
    a = ap.parse_args()
    sel = [int(s) for s in a.only.split(",") if s.strip()] or [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13] + ([9] if a.templates else [])
    T0 = time.time()
    for s in sel:
        t = time.time(); SECTIONS[s](); print(f"   ({time.time()-t:.1f} s)\n")
    nf = sum(1 for r in RESULTS if not r[2])
    print(f"{len(RESULTS)-nf}/{len(RESULTS)} checks passed in {time.time()-T0:.0f} s"
          + ("" if nf == 0 else "  -- FAILED: " + "; ".join(f"{r[0]}:{r[1]}" for r in RESULTS if not r[2])))
    sys.exit(1 if nf else 0)
