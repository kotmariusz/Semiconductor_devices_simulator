"""SemiSim 3D - triangular-prism meshes for the box-method core.

The device cross-section (the plane normal to the prism axis) is covered by a
Delaunay triangulation and the triangles are extruded along the axis over the
tensor mesher's 1-D grid.  The control volumes are the Voronoi cells of the
triangulation times the 1-D cells, so every edge couples its two nodes through
(dual face area / edge length) exactly as on the tensor mesh, and the core runs
the same Scharfetter-Gummel / Poisson discretisation on the prism graph.

Point placement
  * Axis-aligned geometry (box regions, contacts, sheets).  The candidate
    points are the tensor mesher's lines (they already carry every pinned
    interface pair, inversion-layer row and junction refinement).  Each line
    gets a nesting level - level L keeps lines at least g0*2^L apart - and a
    point is kept only where the local target size admits both of its lines.
    The target size grows with the distance to the nearest feature (material
    interface, junction, gate, contact edge, sheet), so the fine lines that
    the tensor mesh drags across the whole device survive only near what
    needs them.  Near a feature the points ARE the tensor points: interfaces
    and boundary layers are identical in both mesh types.
  * Shaped regions and box contacts (Shape: convex outlines with rounded
    corners, circles; offset families for conformal gate stacks) get rows of
    points along the outline normals: a node pair straddling every interface
    (the Voronoi face between the two IS the interface), a geometric boundary
    layer on the semiconductor side, the metal row on a gate surface.  All
    members of an offset family share the normals, so the layers of a gate
    stack stay aligned (tunnelling paths run straight through them).  Tensor
    points inside the band are dropped.

Units: nm inside this module, metres at the core interface.
"""
import numpy as np
import matplotlib.tri as mtri

AXES = "xyz"
INPLANE = {0: (1, 2), 1: (0, 2), 2: (0, 1)}      # prism axis -> in-plane (u, v) axes
BL_H0, BL_G = 0.5, 1.3                           # boundary layer: first step (nm), growth
DELTA_SEMI = 0.15                                # node pair half-distance at semiconductor interfaces (nm)


# ══════════════════════════════════════════════════════════════════════
#  SHAPES
# ══════════════════════════════════════════════════════════════════════
class Shape:
    """Outline of a region or of a box contact in the plane normal to `axis`
    ('x', 'y' or 'z'): a convex polygon (corner coordinates in nm, the two
    other axes in x-y-z order - for axis 'x' the corners are (y, z)) with an
    optional fillet radius at every corner.  The owner occupies its box
    INTERSECTED with the outline (the box clips it, e.g. a fin shell at the
    STI surface); with hole=True the owner is its box MINUS the outline (a
    gate electrode wrapped around a fin or a nanowire).

    offset(t) returns the parallel outline at distance t (outwards for t > 0;
    corner radii become r + t, exactly like a conformal film).  Outlines of
    one offset family are meshed together: their rows share the normals, so
    the layers of a gate stack stay aligned node by node."""

    def __init__(self, axis, pts, radii=None, hole=False, base=None, off=0.):
        axis = str(axis).lower()
        if axis not in AXES:
            raise ValueError("shape axis must be 'x', 'y' or 'z'")
        P = np.array(pts, dtype=float).reshape(-1, 2)
        n = len(P)
        if n < 3:
            raise ValueError("a shape needs at least 3 corners")
        R = np.zeros(n) if radii is None else np.array(radii, dtype=float).reshape(n)
        a2 = np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])
        if a2 < 0:                                   # make it counter-clockwise
            P = P[::-1].copy(); R = R[::-1].copy()
        self.axis, self.P, self.R, self.hole = axis, P, np.maximum(R, 0.), bool(hole)
        self.base = base if base is not None else (axis, P.copy(), self.R.copy())
        self.off = float(off)
        self._geom()

    # -- construction helpers
    @staticmethod
    def circle(axis, uc, vc, radius, hole=False):
        """Circle of the given radius (nm) centred at (uc, vc)."""
        r = float(radius)
        return Shape(axis, [(uc - r, vc - r), (uc + r, vc - r), (uc + r, vc + r), (uc - r, vc + r)],
                     [r] * 4, hole)

    @staticmethod
    def rounded_rect(axis, u0, v0, u1, v1, radii=0.):
        """Rectangle [u0,u1] x [v0,v1] with corner radii (scalar or 4 values,
        corners in the order (u0,v0), (u1,v0), (u1,v1), (u0,v1))."""
        rr = [radii] * 4 if np.isscalar(radii) else list(radii)
        return Shape(axis, [(u0, v0), (u1, v0), (u1, v1), (u0, v1)], rr)

    def offset(self, t):
        """Parallel outline at distance t (nm, > 0 outwards)."""
        t = float(t)
        n1, n2 = self.n_in, self.n_out
        P = self.P + t * (n1 + n2) / (1. + np.sum(n1 * n2, axis=1))[:, None]
        return Shape(self.axis, P, np.maximum(self.R + t, 0.), self.hole, self.base, self.off + t)

    def outside(self):
        """The same outline as a hole (owner = its box minus the outline)."""
        return Shape(self.axis, self.P, self.R, True, self.base, self.off)

    def to_dict(self):
        ax, bP, bR = self.base
        return dict(axis=self.axis, pts=self.P.tolist(), radii=self.R.tolist(), hole=self.hole,
                    base=dict(pts=np.asarray(bP).tolist(), radii=np.asarray(bR).tolist()), off=self.off)

    @staticmethod
    def from_dict(d):
        if d is None or isinstance(d, Shape):
            return d
        b = d.get("base")
        base = (d["axis"], np.array(b["pts"], float), np.array(b["radii"], float)) if b else None
        return Shape(d["axis"], d["pts"], d.get("radii"), d.get("hole", False), base, d.get("off", 0.))

    def family_key(self):
        ax, bP, bR = self.base
        return (ax, tuple(np.round(np.asarray(bP), 6).ravel()), tuple(np.round(np.asarray(bR), 6)))

    def describe(self):
        k = "hole " if self.hole else ""
        if np.allclose(self.R, self.R[0]) and len(self.P) == 4 and self.R[0] > 0 and \
                np.allclose(np.ptp(self.P, axis=0), 2 * self.R[0]):
            c = self.P.mean(axis=0)
            return f"{k}circle r={self.R[0]:g} at ({c[0]:g},{c[1]:g}) nm ⊥{self.axis}"
        return f"{k}{len(self.P)}-corner outline ⊥{self.axis}" + (" (rounded)" if np.any(self.R > 0) else "")

    # -- geometry
    def _geom(self):
        P, R = self.P, self.R
        n = len(P)
        D = np.roll(P, -1, axis=0) - P                  # edge k: P[k] -> P[k+1]
        L = np.hypot(D[:, 0], D[:, 1])
        if np.any(L < 1e-9):
            raise ValueError("shape has coincident corners")
        d = D / L[:, None]
        nrm = np.c_[d[:, 1], -d[:, 0]]                  # outward normals (CCW)
        cr = d[:, 0] * np.roll(d, -1, axis=0)[:, 1] - d[:, 1] * np.roll(d, -1, axis=0)[:, 0]
        if np.any(cr < -1e-9):
            raise ValueError("shape outline must be convex")
        self.d_in, self.d_out = np.roll(d, 1, axis=0), d          # edge into / out of corner k
        self.n_in, self.n_out = np.roll(nrm, 1, axis=0), nrm
        cphi = np.clip(np.sum(self.n_in * self.n_out, axis=1), -1., 1.)
        self.phi = np.arccos(cphi)                                 # turning angle at corner k
        self.tau = R * np.tan(0.5 * self.phi)
        self.C = P - self.tau[:, None] * self.d_in - R[:, None] * self.n_in   # fillet centres
        self.T_in = P - self.tau[:, None] * self.d_in
        self.T_out = P + self.tau[:, None] * self.d_out
        free = L - self.tau - np.roll(self.tau, -1)
        if np.any(free < -1e-6 * np.maximum(L, 1.)):
            raise ValueError("fillet radii too large for the shape's edges")
        self.Lfree = np.maximum(free, 0.)

    def pieces(self):
        """Boundary as ('A', centre, r, a0, a1) arcs (angles of the outward
        normal) and ('L', p0, p1, normal) lines, counter-clockwise."""
        out = []
        for k in range(len(self.P)):
            a0 = np.arctan2(self.n_in[k, 1], self.n_in[k, 0])
            out.append(("A", self.C[k], self.R[k], a0, a0 + self.phi[k]))
            out.append(("L", self.T_out[k], self.T_in[(k + 1) % len(self.P)], self.n_out[k]))
        return out

    def polyline(self, sag=2e-4):
        """Dense closed polygon on the true outline (sagitta <= sag nm)."""
        pts = []
        for pc in self.pieces():
            if pc[0] == "A":
                _, c, r, a0, a1 = pc
                if r <= 0 or a1 - a0 <= 0:
                    pts.append(c + r * np.array([np.cos(a0), np.sin(a0)])); continue
                step = 2 * np.arccos(max(1. - sag / r, -1.)) if r > sag else a1 - a0
                m = max(1, int(np.ceil((a1 - a0) / max(step, 1e-6))))
                a = a0 + (a1 - a0) * np.arange(m) / m
                pts.extend(c + r * np.c_[np.cos(a), np.sin(a)])
            else:
                pts.append(pc[1])
        return np.array(pts)

    def contains(self, U, V, tol=1e-6):
        """Points (nm arrays) inside the outline, boundary within tol included
        (tol < 0: strictly inside by |tol|)."""
        U = np.asarray(U, float); V = np.asarray(V, float)
        ok = np.ones(np.broadcast(U, V).shape, bool)
        for k in range(len(self.P)):                        # half-planes of the sharp polygon
            nx, ny = self.n_out[k]
            ok &= (U - self.P[k, 0]) * nx + (V - self.P[k, 1]) * ny <= tol
        for k in range(len(self.P)):                        # rounded corners cut away
            if self.R[k] <= 0:
                continue
            cx, cy = self.C[k]
            du, dv = U - cx, V - cy
            wedge = (du * self.n_in[k, 0] + dv * self.n_in[k, 1] > 0) & \
                    (du * self.n_out[k, 0] + dv * self.n_out[k, 1] > 0)
            ok &= ~(wedge & (np.hypot(du, dv) > self.R[k] + tol))
        return ok

    def owner_mask(self, U, V, tol=1e-6):
        """Where the OWNER of this shape is (inside, or outside for a hole;
        a hole's owner includes the outline itself: a gate's metal surface)."""
        return ~self.contains(U, V, -tol) if self.hole else self.contains(U, V, tol)

    def bbox(self):
        pl = self.polyline(1e-3)
        return pl.min(axis=0), pl.max(axis=0)


def _root_shape(sh):
    ax, P, R = sh.base
    return Shape(ax, P, R)


# ══════════════════════════════════════════════════════════════════════
#  REGIONS / CONTACTS AS POINT FUNCTIONS
# ══════════════════════════════════════════════════════════════════════
def _halfopen(X, lo, hi, L, tol):
    """lo <= X < hi (the tensor mesher's rule); an interval reaching the far
    domain edge includes it."""
    m = X >= lo - tol
    return m & (X <= hi + tol) if hi >= L - tol else m & (X < hi - tol)


def _shape_mask(sh, X, Y, Z, tol):
    ai = AXES.index(sh.axis)
    iu, iv = INPLANE[ai]
    C = (X, Y, Z)
    return sh.owner_mask(C[iu], C[iv], tol)


def region_index(regs, X, Y, Z, L, tol=None):
    """Index of the region owning each point (the last one containing it,
    half-open boxes, intersected with the region's shape); -1 = none."""
    X = np.asarray(X, float); Y = np.asarray(Y, float); Z = np.asarray(Z, float)
    tol = 1e-7 * max(L) if tol is None else tol
    own = np.full(np.broadcast(X, Y, Z).shape, -1, np.int32)
    for k, r in enumerate(regs):
        m = _halfopen(X, r.x0, r.x1, L[0], tol) & _halfopen(Y, r.y0, r.y1, L[1], tol) & \
            _halfopen(Z, r.z0, r.z1, L[2], tol)
        sh = getattr(r, "shape", None)
        if sh is not None and m.any():
            m &= _shape_mask(sh, X, Y, Z, 1e-6)
        own[m] = k
    return own


def con_box_nm(c, L):
    """(x0,x1,y0,y1,z0,z1) nm of a box contact."""
    return (c.i0p * L[0], c.i1p * L[0], c.j0p * L[1], c.j1p * L[1], c.k0p * L[2], c.k1p * L[2])


def box_contact_mask(c, X, Y, Z, L, tol=None):
    """Points inside box contact c - half-open box as on the tensor mesh (the
    node ON the far face belongs to what follows), cut to its shape (a shape
    includes its own outline: the metal surface of a wrapped gate)."""
    tol = 1e-7 * max(L) if tol is None else tol
    x0, x1, y0, y1, z0, z1 = con_box_nm(c, L)
    m = _halfopen(X, x0, x1, L[0], tol) & _halfopen(Y, y0, y1, L[1], tol) & _halfopen(Z, z0, z1, L[2], tol)
    sh = getattr(c, "shape", None)
    if sh is not None and m.any():
        m &= _shape_mask(sh, X, Y, Z, 1e-6)
    return m


def face_axes(face):
    """(normal axis, first in-plane axis, second in-plane axis) of a face
    contact (the i and j ranges of Con)."""
    return {0: (0, 1, 2), 1: (0, 1, 2), 2: (1, 0, 2), 3: (1, 0, 2), 4: (2, 0, 1), 5: (2, 0, 1)}[face]


def face_range_nm(c, L):
    """Face contact: (normal axis, side, (a-axis, a0, a1), (b-axis, b0, b1))."""
    na, a_ax, b_ax = face_axes(c.face)
    return (na, c.face % 2, (a_ax, c.i0p * L[a_ax], c.i1p * L[a_ax]), (b_ax, c.j0p * L[b_ax], c.j1p * L[b_ax]))


# ══════════════════════════════════════════════════════════════════════
#  1-D LINE LEVELS AND SMALL HELPERS
# ══════════════════════════════════════════════════════════════════════
def line_levels(c):
    """Nesting levels of the sorted 1-D grid c: level L keeps lines at least
    g0*2^L apart (greedy from the low end, domain ends always kept);
    lev[i] = the coarsest level that still contains line i."""
    c = np.asarray(c, float)
    n = len(c)
    g0 = max(float(np.min(np.diff(c))), 1e-9)
    lev = np.zeros(n, np.int64)
    alive = np.arange(n)
    L = 0
    while len(alive) > 2 and L < 60:
        L += 1
        G = g0 * 2.0 ** L
        keep = [alive[0]]
        for i in alive[1:-1]:
            if c[i] - c[keep[-1]] >= G and c[-1] - c[i] >= 0.5 * G:
                keep.append(i)
        keep.append(alive[-1])
        alive = np.array(keep)
        lev[alive] = L
    lev[0] = lev[-1] = 10 ** 6
    return lev, g0


def _seg_dist(P, a, b):
    """Distance of points P (n,2) to segment a-b."""
    ab = b - a
    L2 = float(ab @ ab)
    if L2 <= 0:
        return np.hypot(P[:, 0] - a[0], P[:, 1] - a[1])
    t = np.clip(((P[:, 0] - a[0]) * ab[0] + (P[:, 1] - a[1]) * ab[1]) / L2, 0., 1.)
    return np.hypot(P[:, 0] - a[0] - t * ab[0], P[:, 1] - a[1] - t * ab[1])


class _Buckets:
    """Uniform-grid bucketing of 2-D points for radius queries."""

    def __init__(self, Q, cell):
        self.Q = np.asarray(Q, float)
        self.cell = float(cell)
        self.map = {}
        if len(self.Q):
            key = np.floor(self.Q / self.cell).astype(np.int64)
            order = np.lexsort((key[:, 1], key[:, 0]))
            ks = key[order]
            brk = np.nonzero(np.any(np.diff(ks, axis=0) != 0, axis=1))[0] + 1
            for grp in np.split(order, brk):
                self.map[(int(key[grp[0], 0]), int(key[grp[0], 1]))] = grp

    def nearest(self, P, reach):
        """Distance to (and index of) the nearest stored point, searched in
        growing square windows up to `reach`; inf / -1 beyond it."""
        P = np.asarray(P, float)
        dist = np.full(len(P), np.inf); idx = np.full(len(P), -1, np.int64)
        if not self.map or not len(P):
            return dist, idx
        key = np.floor(P / self.cell).astype(np.int64)
        todo = np.arange(len(P))
        m = 1
        while len(todo):
            order = todo[np.lexsort((key[todo, 1], key[todo, 0]))]
            ks = key[order]
            brk = np.nonzero(np.any(np.diff(ks, axis=0) != 0, axis=1))[0] + 1
            for grp in np.split(order, brk):
                ki, kj = int(key[grp[0], 0]), int(key[grp[0], 1])
                cand = [self.map[(ki + a, kj + b)] for a in range(-m, m + 1) for b in range(-m, m + 1)
                        if (ki + a, kj + b) in self.map]
                if not cand:
                    continue
                cand = np.concatenate(cand)
                for g0 in range(0, len(grp), 512):            # bounded distance blocks
                    gg = grp[g0:g0 + 512]
                    dd = np.hypot(P[gg, None, 0] - self.Q[None, cand, 0], P[gg, None, 1] - self.Q[None, cand, 1])
                    j = np.argmin(dd, axis=1)
                    dist[gg] = dd[np.arange(len(gg)), j]; idx[gg] = cand[j]
            # exact once the nearest lies within the window's inner radius m*cell
            todo = todo[~(dist[todo] <= m * self.cell)]
            if (m - 1) * self.cell > reach:
                break
            m *= 2
        far = dist > reach
        dist[far] = np.inf; idx[far] = -1
        return dist, idx


def _dedupe(P, tol):
    key = np.round(P / tol).astype(np.int64)
    _, first = np.unique(key, axis=0, return_index=True)
    return P[np.sort(first)]


# ══════════════════════════════════════════════════════════════════════
#  THE PRISM MESH
# ══════════════════════════════════════════════════════════════════════
class PrismMesh:
    """Triangular-prism mesh of a device dict t (template format).
    lines: (xs, ys, zs) in metres from the tensor mesher (main.make_grid).
    insulators: set of insulator material names.  g_feat / g_junc: growth
    of the target spacing with the distance to an interface / a junction;
    thin=False keeps every tensor point (box geometry then gives the tensor
    mesh itself, triangulated - used to validate the prism path)."""

    def __init__(self, t, lines, axis, insulators, g_feat=0.3, g_junc=0.2, thin=True):
        if not thin:                        # keep every tensor point (validation: = tensor mesh)
            g_feat = g_junc = 0.
        self.t = t
        self.L = (float(t["Lx"]), float(t["Ly"]), float(t["Lz"]))
        self.ax = AXES.index(axis) if isinstance(axis, str) else int(axis)
        self.iu, self.iv = INPLANE[self.ax]
        self.ins = set(insulators)
        C = [np.asarray(c, float) * 1e9 for c in lines]
        self.lines = C
        self.U, self.V, self.A = C[self.iu], C[self.iv], C[self.ax]
        self.Lu, self.Lv, self.La = self.L[self.iu], self.L[self.iv], self.L[self.ax]
        self.regs, self.cons = t["regs"], t["cons"]
        self.sheets = t.get("sheets", []) or []
        self.g_feat, self.g_junc = g_feat, g_junc
        self.tol = 1e-7 * max(self.L)
        self.warn = []
        self._layer_reps()
        self._features()
        self._families()
        self._points()
        self._geometry()
        self._extrude()

    # -- representative layer coordinates (one per interval between region /
    #    contact limits along the prism axis)
    def _layer_reps(self):
        ax = self.ax
        cuts = {0., self.La}
        for r in self.regs:
            cuts.update(((r.x0, r.y0, r.z0)[ax], (r.x1, r.y1, r.z1)[ax]))
        for c in self.cons:
            if c.face == 6:
                b = con_box_nm(c, self.L); cuts.update((b[2 * ax], b[2 * ax + 1]))
        cuts = np.array(sorted(v for v in cuts if -1e-9 <= v <= self.La + 1e-9))
        mids = 0.5 * (cuts[1:] + cuts[:-1])
        mids = mids[np.diff(cuts) > 1e-9]
        self.a_reps = mids if len(mids) else np.array([0.5 * self.La])

    def _xyz(self, u, v, a):
        """In-plane (u, v) and axial a -> X, Y, Z arrays (broadcast)."""
        out = [None, None, None]
        out[self.iu], out[self.iv], out[self.ax] = u, v, a
        return np.broadcast_arrays(*out)

    def identity(self, u, v):
        """Per representative layer: (material id, doping signature, contact id)
        at in-plane points -> arrays (n_layers, n_points)."""
        u = np.asarray(u, float); v = np.asarray(v, float)
        nA = len(self.a_reps)
        uu = np.broadcast_to(u, (nA,) + u.shape); vv = np.broadcast_to(v, (nA,) + u.shape)
        aa = np.broadcast_to(self.a_reps.reshape((nA,) + (1,) * u.ndim), (nA,) + u.shape)
        X, Y, Z = self._xyz(uu, vv, aa)
        own = region_index(self.regs, X, Y, Z, self.L, self.tol)
        mats = np.array([r.mat for r in self.regs] + ["<none>"])
        ins = np.array([r.mat in self.ins for r in self.regs] + [True])
        dop = np.array([0. if (r.mat in self.ins or r.doping <= 1e11) else
                        (r.doping if r.dtype == "n" else -r.doping) for r in self.regs] + [0.])
        mat = mats[own]; semi = ~ins[own]; dsig = dop[own]
        con = np.full(own.shape, -1, np.int32)
        for k, c in enumerate(self.cons):
            if c.face == 6:
                con[box_contact_mask(c, X, Y, Z, self.L, self.tol)] = k
        return mat, semi, dsig, con

    # -- axis-aligned features from region / contact boxes
    def _features(self):
        Lu, Lv = self.Lu, self.Lv
        iu, iv, ax = self.iu, self.iv, self.ax
        tol = self.tol
        segs = []                                  # (orient, c, a, b): 'u' line u=c, v in [a,b]
        rects = []
        for r in self.regs:
            lo = (r.x0, r.y0, r.z0); hi = (r.x1, r.y1, r.z1)
            rects.append((lo[iu], hi[iu], lo[iv], hi[iv]))
        for c in self.cons:
            if c.face == 6:
                b = con_box_nm(c, self.L)
                rects.append((b[2 * iu], b[2 * iu + 1], b[2 * iv], b[2 * iv + 1]))
        for (u0, u1, v0, v1) in rects:
            u0, u1 = max(u0, 0.), min(u1, Lu); v0, v1 = max(v0, 0.), min(v1, Lv)
            if u1 - u0 <= tol or v1 - v0 <= tol:
                continue
            for cc in (u0, u1):
                if tol < cc < Lu - tol: segs.append(("u", cc, v0, v1))
            for cc in (v0, v1):
                if tol < cc < Lv - tol: segs.append(("v", cc, u0, u1))
        # split points: every other line crossing, and crossings with shape outlines
        ucut = sorted({s[1] for s in segs if s[0] == "u"} | {0., Lu})
        vcut = sorted({s[1] for s in segs if s[0] == "v"} | {0., Lv})
        polys = []
        for obj in list(self.regs) + list(self.cons):
            sh = getattr(obj, "shape", None)
            if sh is not None and AXES.index(sh.axis) == ax:
                polys.append(sh.polyline(1e-3))
        pieces = []
        for (o, cc, a, b) in segs:
            cuts = set((vcut if o == "u" else ucut))
            for pl in polys:                         # outline crossings of this line
                q = pl[:, 0] if o == "u" else pl[:, 1]
                w = pl[:, 1] if o == "u" else pl[:, 0]
                q2, w2 = np.roll(q, -1), np.roll(w, -1)
                hit = (q - cc) * (q2 - cc) <= 0
                for i in np.nonzero(hit & (q != q2))[0]:
                    cuts.add(w[i] + (cc - q[i]) * (w2[i] - w[i]) / (q2[i] - q[i]))
            cs = np.array(sorted(x for x in cuts if a + tol < x < b - tol))
            ends = np.concatenate(([a], cs, [b]))
            for p0, p1 in zip(ends[:-1], ends[1:]):
                if p1 - p0 > tol:
                    pieces.append((o, cc, p0, p1))
        pieces = sorted(set(pieces))
        # classify each piece by the material / doping / contact on its two sides
        feats = []
        if pieces:
            eps = max(1e-4, 10 * tol)
            mid = np.array([0.5 * (p[2] + p[3]) for p in pieces])
            cc = np.array([p[1] for p in pieces])
            isu = np.array([p[0] == "u" for p in pieces])
            uA = np.where(isu, cc - eps, mid); vA = np.where(isu, mid, cc - eps)
            uB = np.where(isu, cc + eps, mid); vB = np.where(isu, mid, cc + eps)
            mA, sA, dA, cA = self.identity(uA, vA)
            mB, sB, dB, cB = self.identity(uB, vB)
            difm = np.any(mA != mB, axis=0); difc = np.any(cA != cB, axis=0)
            difd = np.any(dA != dB, axis=0)
            semA = np.any(sA & ~sB & (mA != mB), axis=0); semB = np.any(sB & ~sA & (mA != mB), axis=0)
            for k, p in enumerate(pieces):
                if difm[k] or difc[k]:
                    kind = "if"
                elif difd[k]:
                    kind = "jn"
                else:
                    continue
                o, c0, a0, b0 = p
                if o == "u": p0, p1 = np.array([c0, a0]), np.array([c0, b0])
                else: p0, p1 = np.array([a0, c0]), np.array([b0, c0])
                feats.append((kind, p0, p1))
        # contacts
        for c in self.cons:
            if c.face == 6:
                continue
            na, side, (a_ax, a0, a1), (b_ax, b0, b1) = face_range_nm(c, self.L)
            if na == ax:                               # contact on an end face: its rectangle
                rr = {a_ax: (a0, a1), b_ax: (b0, b1)}
                (u0, u1), (v0, v1) = rr[iu], rr[iv]
                for P0, P1 in (((u0, v0), (u1, v0)), ((u1, v0), (u1, v1)), ((u1, v1), (u0, v1)), ((u0, v1), (u0, v0))):
                    feats.append(("ce", np.array(P0, float), np.array(P1, float)))
                continue
            fc = 0. if side == 0 else self.L[na]
            tang = a_ax if a_ax != ax else b_ax
            t0, t1 = (a0, a1) if tang == a_ax else (b0, b1)
            if na == iu: p0, p1 = np.array([fc, t0]), np.array([fc, t1])
            else: p0, p1 = np.array([t0, fc]), np.array([t1, fc])
            if c.bc in (1, 2) or getattr(c, "metal", ""):
                # gates, Schottky and metal contacts: fine all along the face
                # (inversion layer / depletion barrier under the whole contact)
                feats.append(("gate", p0, p1))
            else:
                feats.append(("ce", p0, p0.copy())); feats.append(("ce", p1, p1.copy()))
        # mesh hints (defects that change the material or carry lines/planes):
        # the outline of their bounding box in the cross-section, refined like
        # an interface
        for lo, hi in (self.t.get("mesh_hints") or []):
            u0, u1, v0, v1 = lo[iu], hi[iu], lo[iv], hi[iv]
            for P0, P1 in (((u0, v0), (u1, v0)), ((u1, v0), (u1, v1)), ((u1, v1), (u0, v1)), ((u0, v1), (u0, v0))):
                feats.append(("if", np.array(P0, float), np.array(P1, float)))
        for s in self.sheets:
            sa = AXES.index(s.axis)
            if sa == ax:
                continue
            oth = [q for q in range(3) if q != sa]
            rng = {oth[0]: (s.a0, s.a1), oth[1]: (s.b0, s.b1)}
            t0, t1 = rng[iv] if sa == iu else rng[iu]
            if sa == iu: p0, p1 = np.array([s.pos, t0]), np.array([s.pos, t1])
            else: p0, p1 = np.array([t0, s.pos]), np.array([t1, s.pos])
            feats.append(("sheet", p0, p1))
        self.feats = feats

    # -- offset families of shapes (curved bands)
    def _families(self):
        fam = {}
        for obj in list(self.regs) + list(self.cons):
            sh = getattr(obj, "shape", None)
            if sh is None:
                continue
            if AXES.index(sh.axis) != self.ax:
                self.warn.append(f"shape of '{getattr(obj, 'label', '?')}' is not normal to the prism "
                                 f"axis - meshed as a staircase")
                continue
            key = sh.family_key()
            f = fam.setdefault(key, dict(root=_root_shape(sh), offs=set()))
            f["offs"].add(round(sh.off, 9))
        self.families = list(fam.values())
        self.band_pts = np.zeros((0, 2))
        self.band_feat = np.zeros((0, 2))
        if not self.families:
            return
        hbase = min(self.Lu / max(len(self.U) - 1, 1), self.Lv / max(len(self.V) - 1, 1))
        allp, featp = [], []
        for f in self.families:
            pts, fp = self._family_rows(f["root"], sorted(f["offs"]), hbase)
            allp.append(pts); featp.append(fp)
        P = np.concatenate(allp) if allp else np.zeros((0, 2))
        inside = (P[:, 0] > -1e-9) & (P[:, 0] < self.Lu + 1e-9) & (P[:, 1] > -1e-9) & (P[:, 1] < self.Lv + 1e-9)
        self.band_pts = _dedupe(P[inside], 1e-5) if inside.any() else np.zeros((0, 2))
        self.band_feat = np.concatenate(featp) if featp else np.zeros((0, 2))

    def _family_rows(self, root, offs, hbase):
        """Rows of points along the root outline's normals for one family."""
        per = sum((pc[2] * (pc[4] - pc[3]) if pc[0] == "A" else np.hypot(*(pc[2] - pc[1])))
                  for pc in root.pieces())
        # tangential spacing: ~48 samples round the outline (fillets at most
        # 15 degrees per step); the boundary layer grows until its steps
        # reach this spacing
        s = float(np.clip(per / 48., 0.4, max(2. * hbase, 0.4)))
        omax = max(max(offs), 0.)
        # samples along the root outline: points, normals, local fillet radius
        P, N, Rl = [], [], []
        for pc in root.pieces():
            if pc[0] == "A":
                _, c, r, a0, a1 = pc
                if a1 - a0 <= 1e-12:
                    continue
                m = max(1, int(np.ceil((a1 - a0) * (r + omax) / s)), int(np.ceil((a1 - a0) / np.radians(15.))))
                a = a0 + (a1 - a0) * np.arange(m) / m
                nn = np.c_[np.cos(a), np.sin(a)]
                P.append(c + r * nn); N.append(nn); Rl.append(np.full(m, r))
            else:
                _, p0, p1, nrm = pc
                ln = np.hypot(*(p1 - p0))
                if ln <= 1e-9:
                    continue
                m = max(1, int(np.ceil(ln / s)))
                w = np.arange(m) / m
                P.append(p0 + w[:, None] * (p1 - p0)); N.append(np.tile(nrm, (m, 1))); Rl.append(np.full(m, np.inf))
        P, N, Rl = np.concatenate(P), np.concatenate(N), np.concatenate(Rl)
        M = len(P)
        offs = np.array(offs, float)
        eps = 1e-4
        # identity on both sides of every offset outline, every layer
        kinds = np.zeros((M, len(offs)), np.int8)        # 0 none, 1 junction, 2 interface, 3 metal face
        blside = np.zeros((M, len(offs), 2), bool)        # boundary layer inward / outward
        metal_out = np.zeros((M, len(offs)), bool)
        for j, o in enumerate(offs):
            qa = P + (o - eps) * N; qb = P + (o + eps) * N
            mA, sA, dA, cA = self.identity(qa[:, 0], qa[:, 1])
            mB, sB, dB, cB = self.identity(qb[:, 0], qb[:, 1])
            difm = np.any(mA != mB, axis=0); difc = np.any(cA != cB, axis=0); difd = np.any(dA != dB, axis=0)
            kinds[:, j] = np.where(difc, 3, np.where(difm, 2, np.where(difd, 1, 0)))
            metal_out[:, j] = np.any((cB >= 0) & (cA != cB), axis=0)
            # semiconductor side of a material interface gets the boundary layer
            hetero = np.any(sA & sB & (mA != mB), axis=0)
            blside[:, j, 0] = np.any(sA & ~sB & (mA != mB), axis=0) | hetero | np.any(sA & (cB >= 0) & (cA < 0), axis=0)
            blside[:, j, 1] = np.any(sB & ~sA & (mA != mB), axis=0) | hetero | np.any(sB & (cA >= 0) & (cB < 0), axis=0)
        # thickness to the opposite side (inward BL limit) - convex outline
        def chord_in(p, n):
            best = np.full(len(p), np.inf)
            for k in range(len(root.P)):
                nk = root.n_out[k]
                den = -(n @ nk)
                num = (root.P[k] - p) @ nk
                t = np.where(den > 1e-12, num / np.where(den > 1e-12, den, 1.), np.inf)
                best = np.minimum(best, np.where(t > 1e-9, t, np.inf))
            return best
        thick = chord_in(P, N)
        pts, prow, featp = [], [], []
        for i in range(M):
            act = np.nonzero(kinds[i] > 0)[0]
            if len(act) == 0:
                continue
            rows = []
            for j in act:
                o = offs[j]
                gaps = [abs(o - offs[q]) for q in act if q != j] or [np.inf]
                # half-distance of the pair: 0.15 nm next to a semiconductor boundary layer
                # (its first cell holds the inversion / accumulation charge), else 0.25 nm
                dl = min(DELTA_SEMI if blside[i, j].any() else 0.25, 0.25 * min(gaps))
                if kinds[i, j] == 3:                          # metal face: metal row on the surface
                    rows += [o, o - 2 * dl] if metal_out[i, j] else [o, o + 2 * dl]
                else:
                    rows += [o - dl, o + dl]
                featp.append(P[i] + o * N[i])
                for sd in (0, 1):
                    if not blside[i, j, sd]:
                        continue
                    sg = -1. if sd == 0 else 1.
                    nxt = [offs[q] for q in act if (offs[q] - o) * sg > 1e-9]
                    lim = (abs(min(nxt, key=lambda q: abs(q - o)) - o) * 0.5) if nxt else np.inf
                    if sg < 0:
                        # inward: at most 45 % of the local thickness, and half a fillet radius
                        # (deeper rows of a fillet converge on its centre: fans of slivers)
                        lim = min(lim, 0.45 * (thick[i] + o), 0.5 * Rl[i] + o if np.isfinite(Rl[i]) else np.inf)
                    lim = min(lim, 15.)
                    d, h = dl, BL_H0
                    while True:
                        d += h
                        if d > lim or h > 1.2 * s:
                            break
                        rows.append(o + sg * d)
                        h *= BL_G
            rows = np.unique(np.round(np.array(rows), 9))
            if np.isfinite(Rl[i]):
                rows = rows[Rl[i] + rows > max(0.25 * s, 0.45 * Rl[i])]   # inner rows of a fillet: not near its centre
            pts.append(P[i][None, :] + rows[:, None] * N[i][None, :]); prow.append(rows)
        if not pts:
            return np.zeros((0, 2)), np.zeros((0, 2))
        Q = np.concatenate(pts); O = np.concatenate(prow)
        # a point of an inward row must still lie |o| inside the outline: near a
        # sharp or tight corner the row of one side runs into the other side
        inner = O < -1e-9
        if inner.any():
            pl = root.polyline(1e-3)
            qi = Q[inner]
            dmin = np.full(len(qi), np.inf)
            for k in range(len(pl)):
                dmin = np.minimum(dmin, _seg_dist(qi, pl[k], pl[(k + 1) % len(pl)]))
            bad = dmin < 0.9 * np.abs(O[inner]) - 1e-6
            keep = np.ones(len(Q), bool); keep[np.nonzero(inner)[0][bad]] = False
            Q = Q[keep]
        self._band_s = getattr(self, "_band_s", []) + [s]
        return Q, np.array(featp)

    # -- point selection
    def _points(self):
        U, V = self.U, self.V
        NU, NV = len(U), len(V)
        levU, g0u = line_levels(U)
        levV, g0v = line_levels(V)
        Ug, Vg = np.meshgrid(U, V)                          # (NV, NU)
        P = np.c_[Ug.ravel(), Vg.ravel()]
        H = np.full(len(P), np.inf)
        for kind, p0, p1 in self.feats:
            g = self.g_junc if kind == "jn" else self.g_feat
            H = np.minimum(H, g * _seg_dist(P, p0, p1))
        hcu = self.Lu / max(self.t["Nx" if self.iu == 0 else ("Ny" if self.iu == 1 else "Nz")] - 1, 1)
        hcv = self.Lv / max(self.t["Nx" if self.iv == 0 else ("Ny" if self.iv == 1 else "Nz")] - 1, 1)
        if len(self.band_feat):
            reach = max(hcu, hcv) / self.g_feat               # farther away the cap rules anyway
            B = _Buckets(self.band_feat, max(0.5, reach / 8.))
            dmin, _ = B.nearest(P, reach)
            H = np.minimum(H, self.g_feat * dmin)
        Hu = np.minimum(H, hcu); Hv = np.minimum(H, hcv)
        Lu_ = np.floor(np.log2(np.maximum(Hu, 1e-30) / g0u)).clip(0, None)
        Lv_ = np.floor(np.log2(np.maximum(Hv, 1e-30) / g0v)).clip(0, None)
        iu_ = np.tile(np.arange(NU), NV); iv_ = np.repeat(np.arange(NV), NU)
        keep = (levU[iu_] >= Lu_) & (levV[iv_] >= Lv_)
        Pk = P[keep]
        self.n_tensor = len(P)
        if len(self.band_pts):
            # drop tensor points in / next to the shape bands
            gu = np.diff(U); gv = np.diff(V)
            hu = np.maximum(np.r_[gu[:1], gu], np.r_[gu, gu[-1:]])
            hv = np.maximum(np.r_[gv[:1], gv], np.r_[gv, gv[-1:]])
            hloc = np.maximum(hu[iu_[keep]], hv[iv_[keep]])
            sband = max(getattr(self, "_band_s", [1.]))
            B = _Buckets(self.band_pts, max(sband, 0.5))
            clear = 0.75 * np.maximum(sband, np.minimum(hloc, 4 * sband))
            dmin, _ = B.nearest(Pk, float(clear.max()) * 1.01)
            # points on the domain edge only give way to a band that nearly touches them
            # (a long bare edge next to a band row makes obtuse boundary triangles)
            onb = (Pk[:, 0] <= 1e-9) | (Pk[:, 0] >= self.Lu - 1e-9) | (Pk[:, 1] <= 1e-9) | (Pk[:, 1] >= self.Lv - 1e-9)
            clear = np.where(onb, 0.3 * clear, clear)
            Pk = Pk[dmin >= clear]
            # band points on the domain edge only if exactly on it
            Pk = np.concatenate([Pk, self.band_pts])
        self.P2 = _dedupe(Pk, 1e-6)

    # -- Delaunay and box-method geometry of the cross-section
    def _geometry(self):
        P = self.P2
        tri = mtri.Triangulation(P[:, 0], P[:, 1])
        T = tri.triangles
        A, B, Cc = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]

        def cot(p, q, r):
            u = q - p; v = r - p
            cr = u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]
            return (u[:, 0] * v[:, 0] + u[:, 1] * v[:, 1]) / np.where(np.abs(cr) > 0, np.abs(cr), 1e-300), np.abs(cr)

        ca, cr = cot(A, B, Cc); cb, _ = cot(B, Cc, A); cc, _ = cot(Cc, A, B)
        area = 0.5 * cr
        good = area > 1e-12 * self.Lu * self.Lv
        if not good.all():
            T = T[good]; A, B, Cc = A[good], B[good], Cc[good]; ca, cb, cc, area = ca[good], cb[good], cc[good], area[good]
        used = np.zeros(len(P), bool); used[T.ravel()] = True
        if not used.all():                      # points dropped by the triangulation (should not happen)
            P = P[used]; remap = -np.ones(len(used), np.int64); remap[used] = np.arange(used.sum())
            T = remap[T]; self.P2 = P
            self.warn.append(f"{(~used).sum()} coincident points removed")
        N2 = len(P)
        ei = np.concatenate([T[:, 1], T[:, 2], T[:, 0]]); ej = np.concatenate([T[:, 2], T[:, 0], T[:, 1]])
        w = 0.5 * np.concatenate([ca, cb, cc])
        lo = np.minimum(ei, ej).astype(np.int64); hi = np.maximum(ei, ej).astype(np.int64)
        key = lo * N2 + hi
        uk, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
        wsum = np.bincount(inv, weights=w)
        E0 = uk // N2; E1 = uk % N2
        Le = np.hypot(*(P[E1] - P[E0]).T)
        dual = wsum * Le                                   # signed dual length (nm)
        self.n_negw = int(np.sum(wsum < -1e-9)); self.n_edges_all = len(uk)
        bnd = cnt == 1
        # control areas: circumcentric (signed Voronoi) per node; mixed area as a fallback
        Ecell = 0.25 * Le * dual
        A2 = np.bincount(E0, Ecell, N2) + np.bincount(E1, Ecell, N2)
        la2 = np.sum((B - Cc) ** 2, 1); lb2 = np.sum((Cc - A) ** 2, 1); lc2 = np.sum((A - B) ** 2, 1)
        obt = (ca < 0) | (cb < 0) | (cc < 0)
        aA = np.where(obt, np.where(ca < 0, area / 2, area / 4), (lb2 * cb + lc2 * cc) / 8)
        aB = np.where(obt, np.where(cb < 0, area / 2, area / 4), (lc2 * cc + la2 * ca) / 8)
        aC = np.where(obt, np.where(cc < 0, area / 2, area / 4), (la2 * ca + lb2 * cb) / 8)
        Amix = np.bincount(T[:, 0], aA, N2) + np.bincount(T[:, 1], aB, N2) + np.bincount(T[:, 2], aC, N2)
        badA = ~(A2 > 0.2 * Amix)
        self.n_area_fix = int(badA.sum())
        A2 = np.where(badA, Amix, A2)
        # boundary dual lengths per domain side (0 u=0, 1 u=Lu, 2 v=0, 3 v=Lv)
        bl = np.zeros((N2, 4))
        tb = 1e-6 * max(self.Lu, self.Lv)
        for s_, (cix, val) in enumerate(((0, 0.), (0, self.Lu), (1, 0.), (1, self.Lv))):
            on = (np.abs(P[E0[bnd], cix] - val) < tb) & (np.abs(P[E1[bnd], cix] - val) < tb)
            e0, e1, le = E0[bnd][on], E1[bnd][on], Le[bnd][on]
            bl[:, s_] += np.bincount(e0, 0.5 * le, N2) + np.bincount(e1, 0.5 * le, N2)
        keepE = dual > 1e-9 * Le
        self.tri = mtri.Triangulation(P[:, 0], P[:, 1], T)
        self.T = T
        self.E0, self.E1, self.Le, self.De = E0[keepE], E1[keepE], Le[keepE], dual[keepE]
        self.A2, self.bl, self.N2 = A2, bl, N2
        # quality
        ang = np.degrees(np.arccos(np.clip(np.concatenate([ca, cb, cc]) / np.sqrt(1 + np.concatenate([ca, cb, cc]) ** 2), -1, 1)))
        self.min_angle = float(np.min(ang)) if len(ang) else 0.
        self.n_obtuse = int(obt.sum())
        self.area_err = float(abs(A2.sum() - self.Lu * self.Lv) / (self.Lu * self.Lv))

    # -- extrusion: layers, nodes, graph
    def _extrude(self):
        A = self.A
        NL = len(A)
        LA = np.empty(NL)
        LA[0] = 0.5 * (A[1] - A[0]); LA[-1] = 0.5 * (A[-1] - A[-2]); LA[1:-1] = 0.5 * (A[2:] - A[:-2])
        self.NL, self.LA = NL, LA
        self.N = self.N2 * NL

    def graph(self):
        """(N, E, rp, cj, h, ar, vol, pos) in SI units for api3_set_mesh_graph."""
        N2, NL = self.N2, self.NL
        E0, E1 = self.E0, self.E1
        ne = len(E0)
        lay = np.arange(NL)
        # in-plane directed edges, every layer
        s_in = (np.concatenate([E0, E1])[None, :] + N2 * lay[:, None]).ravel()
        d_in = (np.concatenate([E1, E0])[None, :] + N2 * lay[:, None]).ravel()
        h_in = np.tile(np.concatenate([self.Le, self.Le]), NL)
        a_in = (np.concatenate([self.De, self.De])[None, :] * self.LA[:, None]).ravel()
        # axial edges
        i2 = np.arange(N2)
        lo = (i2[None, :] + N2 * lay[:-1, None]).ravel(); hi = lo + N2
        dz = np.repeat(np.diff(self.A), N2)
        aa = np.tile(self.A2, NL - 1)
        src = np.concatenate([s_in, lo, hi]); dst = np.concatenate([d_in, hi, lo])
        h = np.concatenate([h_in, dz, dz]) * 1e-9
        ar = np.concatenate([a_in, aa, aa]) * 1e-18
        order = np.lexsort((dst, src))
        src, dst, h, ar = src[order], dst[order], h[order], ar[order]
        N = N2 * NL
        rp = np.zeros(N + 1, np.int32)
        rp[1:] = np.cumsum(np.bincount(src, minlength=N))
        vol = (self.A2[None, :] * self.LA[:, None]).ravel() * 1e-27
        pos = np.zeros((N, 3))
        pos[:, self.iu] = np.tile(self.P2[:, 0], NL)
        pos[:, self.iv] = np.tile(self.P2[:, 1], NL)
        pos[:, self.ax] = np.repeat(self.A, N2)
        return (N, len(src), rp, dst.astype(np.int32), h, ar, vol, pos.ravel() * 1e-9)

    def node_xyz(self):
        """Node coordinates (N,3) in nm, node k = i2 + N2*layer."""
        N2, NL = self.N2, self.NL
        pos = np.zeros((self.N, 3))
        pos[:, self.iu] = np.tile(self.P2[:, 0], NL)
        pos[:, self.iv] = np.tile(self.P2[:, 1], NL)
        pos[:, self.ax] = np.repeat(self.A, N2)
        return pos

    # -- node data
    def node_regions(self):
        X = self.node_xyz()
        return region_index(self.regs, X[:, 0], X[:, 1], X[:, 2], self.L, self.tol)

    def sheet_charge(self):
        """Fixed charge density Nf (m^-3) per node from the sheets."""
        Nf = np.zeros(self.N)
        N2 = self.N2
        for s in self.sheets:
            sa = AXES.index(s.axis)
            sig = s.sigma * 1e4                              # q/m^2
            oth = [q for q in range(3) if q != sa]
            rng = {oth[0]: (s.a0, s.a1), oth[1]: (s.b0, s.b1)}
            if sa == self.ax:
                l = int(np.argmin(np.abs(self.A - s.pos)))
                (u0, u1), (v0, v1) = rng[self.iu], rng[self.iv]
                m = (self.P2[:, 0] >= u0 - self.tol) & (self.P2[:, 0] <= u1 + self.tol) & \
                    (self.P2[:, 1] >= v0 - self.tol) & (self.P2[:, 1] <= v1 + self.tol)
                Nf[np.nonzero(m)[0] + N2 * l] += sig / (self.LA[l] * 1e-9)
                continue
            (a0, a1) = rng[self.ax]
            t0, t1 = rng[self.iv] if sa == self.iu else rng[self.iu]
            t0, t1 = max(t0, 0.), min(t1, self.Lv if sa == self.iu else self.Lu)
            if t1 <= t0:
                continue
            m = int(np.ceil((t1 - t0) / max(0.02, (t1 - t0) / 20000.))) + 1
            tt = t0 + (t1 - t0) * (np.arange(m) + 0.5) / m
            dl = (t1 - t0) / m
            q2 = np.zeros(N2)
            for sg in (-1., 1.):
                off = s.pos + sg * 1e-4
                uu, vv = (np.full(m, off), tt) if sa == self.iu else (tt, np.full(m, off))
                k = self.nearest_node2(uu, vv)
                q2 += np.bincount(k, np.full(m, 0.5 * sig * dl * 1e-9), N2)
            lays = np.nonzero((self.A >= a0 - self.tol) & (self.A <= a1 + self.tol))[0]
            for l in lays:
                Nf[N2 * l:N2 * (l + 1)] += q2 / (self.A2 * 1e-18)
        return Nf

    def nearest_node2(self, u, v):
        """Nearest 2-D node (the Voronoi cell) of in-plane points."""
        u = np.asarray(u, float); v = np.asarray(v, float)
        if not hasattr(self, "_nb"):
            self._nb = _Buckets(self.P2, max(0.05, 2 * np.sqrt(self.Lu * self.Lv / max(self.N2, 1))))
        d, k = self._nb.nearest(np.c_[u, v], 2 * max(self.Lu, self.Lv))
        return k

    # -- contacts
    def contact_nodes(self, c):
        """(node indices, areas m^2) of contact c on the prism mesh."""
        N2, NL = self.N2, self.NL
        tol = self.tol
        if c.face == 6:
            X = self.node_xyz()
            m = box_contact_mask(c, X[:, 0], X[:, 1], X[:, 2], self.L, tol)
            k = np.nonzero(m)[0]
            if len(k) == 0:                    # too thin: the node nearest the box centre
                b = con_box_nm(c, self.L)
                ctr = np.array([0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5])])
                k = np.array([int(np.argmin(np.sum((X - ctr) ** 2, axis=1)))])
            return k.astype(np.int32), np.zeros(len(k))
        na, side, (a_ax, a0, a1), (b_ax, b0, b1) = face_range_nm(c, self.L)

        def span(x, lo, hi, L):
            m = _halfopen(x, lo, hi, L, tol)
            if not m.any():
                m = np.zeros(len(x), bool); m[int(np.argmin(np.abs(x - 0.5 * (lo + hi))))] = True
            return m
        if na == self.ax:
            l = 0 if side == 0 else NL - 1
            rr = {a_ax: (a0, a1), b_ax: (b0, b1)}
            (u0, u1), (v0, v1) = rr[self.iu], rr[self.iv]
            m = span(self.P2[:, 0], u0, u1, self.Lu) & span(self.P2[:, 1], v0, v1, self.Lv)
            i2 = np.nonzero(m)[0]
            return (i2 + N2 * l).astype(np.int32), self.A2[i2] * 1e-18
        cix = 0 if na == self.iu else 1
        fc = 0. if side == 0 else self.L[na]
        sd = 2 * cix + side
        onb = np.abs(self.P2[:, cix] - fc) < 1e-6 * max(self.Lu, self.Lv)
        tang = 1 - cix
        tax = self.iu if tang == 0 else self.iv
        t0, t1 = (a0, a1) if tax == a_ax else (b0, b1)
        l0, l1 = (a0, a1) if self.ax == a_ax else (b0, b1)
        Lt = self.Lu if tang == 0 else self.Lv
        cand = np.nonzero(onb)[0]
        mt = span(self.P2[cand, tang], t0, t1, Lt)
        i2 = cand[mt]
        ls = np.nonzero(span(self.A, l0, l1, self.La))[0]
        k = (i2[None, :] + N2 * ls[:, None]).ravel()
        ar = (self.bl[i2, sd][None, :] * self.LA[ls][:, None]).ravel() * 1e-18
        return k.astype(np.int32), ar

    # -- interface walls (MLDA, Lombardi): two per node
    def walls(self, kind_nodes, robin):
        """kind_nodes: per node 0 semiconductor, 1 insulator, 2 metal (as the
        core will see it).  robin: list of (node array, face normal axis, side)
        for oxide (Robin) gates.  Returns the 6N array for api3_set_walls."""
        N2, NL = self.N2, self.NL
        W = np.ones((self.N, 6))
        P = self.P2
        E0, E1, De = self.E0, self.E1, self.De
        # in-plane neighbours per node for cell extents
        nb_i = np.concatenate([E0, E1]); nb_j = np.concatenate([E1, E0])
        order = np.argsort(nb_i, kind="stable"); nb_i, nb_j = nb_i[order], nb_j[order]
        rp = np.zeros(N2 + 1, np.int64); rp[1:] = np.cumsum(np.bincount(nb_i, minlength=N2))
        rob2 = {}
        for nodes, na, side in robin:
            for k in nodes:
                l, i = divmod(int(k), N2)
                rob2.setdefault(l, []).append((i, na, side))
        cache = {}
        for l in range(NL):
            kd = kind_nodes[N2 * l:N2 * (l + 1)]
            rl = tuple(sorted(rob2.get(l, [])))
            key = (kd.tobytes(), rl)
            if key in cache:
                W[N2 * l:N2 * (l + 1)] = cache[key]
                continue
            Wl = np.ones((N2, 6))
            semi = kd == 0
            # facets: semiconductor-insulator edges (midpoint, normal semi -> insulator, dual length)
            fa = semi[E0] & (kd[E1] == 1); fb = semi[E1] & (kd[E0] == 1)
            fs = np.concatenate([E0[fa], E1[fb]]); fi = np.concatenate([E1[fa], E0[fb]])
            fl = np.concatenate([De[fa], De[fb]])
            fm = 0.5 * (P[fs] + P[fi])
            fn = (P[fi] - P[fs]); fn /= np.hypot(fn[:, 0], fn[:, 1])[:, None]
            # Robin gate faces: facet at the boundary node, outward normal
            if rl:
                ri = np.array([i for i, _, _ in rl])
                rn = np.zeros((len(rl), 2))
                for q, (i, na, side) in enumerate(rl):
                    cix = 0 if na == self.iu else 1
                    rn[q, cix] = 1. if side == 1 else -1.
                rlen = np.array([max(self.bl[i].max(), 1e-3) for i, _, _ in rl])
                fm = np.concatenate([fm, P[ri]]); fn = np.concatenate([fn, rn]); fl = np.concatenate([fl, rlen])
            if len(fm) and semi.any():
                si = np.nonzero(semi)[0]
                reach = 60.
                B = _Buckets(fm, 10.)
                # candidate facets per node via buckets (reach 60 nm)
                best = [[] for _ in range(len(si))]
                m = int(np.ceil(reach / B.cell))
                bkey = np.floor(P[si] / B.cell).astype(np.int64)
                for q, (ki, kj) in enumerate(bkey):
                    cand = [B.map[(ki + a, kj + b)] for a in range(-m, m + 1) for b in range(-m, m + 1)
                            if (ki + a, kj + b) in B.map]
                    if not cand:
                        continue
                    cand = np.concatenate(cand)
                    rel = fm[cand] - P[si[q]]
                    z = np.sum(rel * fn[cand], axis=1)
                    lat = np.hypot(rel[:, 0] - z * fn[cand, 0], rel[:, 1] - z * fn[cand, 1])
                    ok = (z > -1e-9) & (z < reach) & (lat <= 0.5 * fl[cand] * (1 + 1e-6) + 1e-6)
                    if not ok.any():
                        continue
                    c2 = cand[ok]; z2 = z[ok]
                    o = np.argsort(z2)
                    walls = [(z2[o[0]], c2[o[0]])]
                    n1 = fn[c2[o[0]]]
                    for oo in o[1:]:
                        if fn[c2[oo]] @ n1 < 0.7:
                            walls.append((z2[oo], c2[oo])); break
                    k = si[q]
                    for w_, (zz, f) in enumerate(walls):
                        nn = fn[f]
                        # cell extent along the wall normal from the in-plane neighbours
                        nbj = nb_j[rp[k]:rp[k + 1]]
                        dv = P[nbj] - P[k]
                        ln = np.hypot(dv[:, 0], dv[:, 1])
                        cs = (dv @ nn) / np.where(ln > 0, ln, 1.)
                        tw = (dv @ nn)[cs > 0.7]; aw = -(dv @ nn)[cs < -0.7]
                        za = max(zz - 0.5 * tw.min(), 0.) if len(tw) else max(zz, 0.)
                        zb = zz + 0.5 * aw.min() if len(aw) else zz
                        Wl[k, 3 * w_:3 * w_ + 3] = (zz * 1e-9, za * 1e-9, zb * 1e-9)
            cache[key] = Wl
            W[N2 * l:N2 * (l + 1)] = Wl
        return W.ravel()

    # -- display: barycentric resampling onto the tensor lines
    def display_setup(self):
        if hasattr(self, "_disp"):
            return self._disp
        Ug, Vg = np.meshgrid(self.U, self.V)
        pu, pv = Ug.ravel(), Vg.ravel()
        finder = self.tri.get_trifinder()
        ti = finder(pu, pv)
        miss = ti < 0
        if miss.any():                                     # nudge towards the domain centre
            cu, cv = 0.5 * self.Lu, 0.5 * self.Lv
            ti[miss] = finder(pu[miss] + 1e-7 * np.sign(cu - pu[miss]), pv[miss] + 1e-7 * np.sign(cv - pv[miss]))
        nodes = np.zeros((len(pu), 3), np.int64); wts = np.zeros((len(pu), 3))
        ok = ti >= 0
        T = self.T[ti[ok]]
        a, b, c = self.P2[T[:, 0]], self.P2[T[:, 1]], self.P2[T[:, 2]]
        det = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
        l1 = ((pu[ok] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (pv[ok] - a[:, 1])) / det
        l2 = ((b[:, 0] - a[:, 0]) * (pv[ok] - a[:, 1]) - (pu[ok] - a[:, 0]) * (b[:, 1] - a[:, 1])) / det
        W3 = np.clip(np.c_[1 - l1 - l2, l1, l2], 0., 1.)
        W3 /= W3.sum(axis=1)[:, None]
        nodes[ok] = T; wts[ok] = W3
        if (~ok).any():
            k = self.nearest_node2(pu[~ok], pv[~ok])
            nodes[~ok] = k[:, None]; wts[~ok] = [1., 0., 0.]
        dom = nodes[np.arange(len(pu)), np.argmax(wts, axis=1)]
        self._disp = (nodes, wts, dom)
        return self._disp

    def to_display(self, vals, mode="lin", mask=None, mat=None):
        """Node values (N) -> array on the tensor lines, shape (Nz, Ny, Nx).
        mode: 'lin' linear, 'log' log10-linear (positive data), 'dom' the
        dominant node.  mask (N bool): only these nodes contribute (NaN where
        none does).  mat (N int): only nodes of the dominant node's material."""
        nodes, wts, dom = self.display_setup()
        N2, NL = self.N2, self.NL
        Vn = np.asarray(vals, float).reshape(NL, N2)
        if mode == "dom":
            out = Vn[:, dom]
        else:
            X = np.log10(np.maximum(Vn, 1e-300)) if mode == "log" else Vn
            w = np.broadcast_to(wts, (NL,) + wts.shape).copy()
            if mask is not None:
                w *= np.asarray(mask, bool).reshape(NL, N2)[:, nodes]
            if mat is not None:
                Mn = np.asarray(mat).reshape(NL, N2)
                w *= Mn[:, nodes] == Mn[:, dom][:, :, None]
            Xv = X[:, nodes]
            w *= np.isfinite(Xv)                            # NaN node values (e.g. Ec on metal) don't count
            sw = w.sum(axis=2)
            Xv = np.where(w > 0, Xv, 0.)
            out = np.where(sw > 1e-12, np.sum(w * Xv, axis=2) / np.where(sw > 1e-12, sw, 1.), np.nan)
            if mode == "log":
                out = 10.0 ** out
        D = out.reshape(NL, len(self.V), len(self.U))      # (a, v, u)
        if self.ax == 2:
            return D                                        # (z, y, x)
        if self.ax == 0:
            return np.ascontiguousarray(D.transpose(1, 2, 0))   # (v=z, u=y, a=x)
        return np.ascontiguousarray(D.transpose(1, 0, 2))       # (v=z, a=y, u=x)

    # -- summary
    def stats(self):
        return dict(N=self.N, N2=self.N2, NL=self.NL, ntri=len(self.T), n_tensor2=self.n_tensor,
                    min_angle=self.min_angle, obtuse=self.n_obtuse, negw=self.n_negw,
                    area_err=self.area_err, area_fix=self.n_area_fix, axis=AXES[self.ax])
