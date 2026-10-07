"""
SemiSim 3D - defects.

A defect is a physical effect (a preset from the library below, its
parameters editable) placed in a shape:

* inclusions: the shape becomes another material - a void, an oxide
  precipitate or a metal/silicide particle - and its surface carries
  interface states (a metal particle: states that pin the Fermi level at its
  Schottky barrier, so it depletes its surroundings like a floating Schottky
  contact and recombines like one);
* point defects / contamination: deep levels (SRH recombination through the
  level + the charge the level traps) with a volume density, in a shape or
  everywhere in the host material; radiation damage from a fluence;
* extended defects: dislocations (deep levels along a line), stacking faults
  and grain boundaries (deep levels on a plane), micropipes (a hollow core);
* interfaces: interface-trap distributions D_it at semiconductor/insulator
  interfaces (and under face gates) and fixed oxide charge.

The library values are representative literature values; capture cross
sections of deep levels are uncertain by up to ~10x (more for Pt), so every
number can be edited.  Placement: manual, or random scatter (count or
density, a seed) - scattered defects are ordinary defects afterwards.

The module maps defects onto a mesh (tensor or prism graph) as node-level
data for the core: material overrides, deep-level species with per-node
densities, and fixed charge.  It has no GUI and no core dependency.
"""
import math
import numpy as np

AXES = "xyz"
Q = 1.602176634e-19

# --------------------------------------------------------------- library
# Trap level: type 'A' acceptor (- when filled), 'D' donor (+ when empty),
# 'N' neutral recombination centre; ref 'Ec' (E below Ec), 'Ev' (E above
# Ev), 'Ei' (E above midgap); capture cross sections in cm^2; w = weight
# (fraction of the defect density, or introduction rate eta in cm^-1 for
# fluence presets).
class Level:
    __slots__ = ("type", "ref", "E", "sn", "sp", "w")
    def __init__(self, type, ref, E, sn, sp, w=1.0):
        self.type, self.ref, self.E, self.sn, self.sp, self.w = type, ref, float(E), float(sn), float(sp), float(w)
    def key(self):
        return (self.type, self.ref, round(self.E, 6), self.sn, self.sp)
    def text(self):
        side = {"Ec": f"Ec − {self.E:.2f}", "Ev": f"Ev + {self.E:.2f}", "Ei": f"Ei {self.E:+.2f}"}[self.ref]
        nm = {"A": "acceptor", "D": "donor", "N": "neutral"}[self.type]
        return f"{nm} {side} eV, σn {self.sn:.1e}, σp {self.sp:.1e} cm²"

L = Level
# geom: what the density means and how the defect acts
#   inclusion  material change in the shape + surface states (dit / pin)
#   volume     deep levels, density in cm^-3, in the shape
#   uniform    the same, by default everywhere in the host material
#   fluence    radiation damage: level densities = eta x fluence (cm^-2)
#   line       deep levels along the cylinder axis, density in cm^-1
#   plane      deep levels on the slab, density in cm^-2 (or a D_it band)
#   interface  D_it (cm^-2 eV^-1) and fixed charge (cm^-2) at
#              semiconductor/insulator interfaces inside the shape
#   charge     fixed charge, cm^-3, in the shape
LIB = {}
def _p(name, cat, geom, desc, **kw):
    d = dict(cat=cat, geom=geom, desc=desc, levels=[], material=None, density=0., dit=0., dit_sigma=1e-15,
             pin=False, metal="", charge=0., shape="sphere", size=(50., 50., 50.), axis="y", host="")
    d.update(kw); LIB[name] = d

# ---- voids and inclusions
_p("Void", "Voids & inclusions", "inclusion",
   "An empty cavity (ε = 1) with an unpassivated inner surface: interface states D_it ~1e12 cm⁻²eV⁻¹ "
   "(amphoteric, σ 1e-15 cm²) that recombine carriers and pin the surface potential. Blocks current, "
   "concentrates the field at its rim.", material="Void", dit=1e12, shape="sphere", size=(100., 100., 100.))
_p("Oxide precipitate (SiO2)", "Voids & inclusions", "inclusion",
   "SiO2 inclusion (oxygen precipitate in Czochralski Si, 10–200 nm) with its interface states "
   "(D_it ~5e11 cm⁻²eV⁻¹): a known recombination/generation centre in CZ wafers.",
   material="SiO2", dit=5e11, shape="sphere", size=(80., 80., 80.))
_p("Metal particle (silicide)", "Voids & inclusions", "inclusion",
   "A metallic inclusion (silicide precipitate, metal particle): a floating equipotential (ε → ∞, "
   "zero net charge) whose surface states pin the Fermi level at the metal's Schottky barrier "
   "(φB from the metal library, 5e13 cm⁻² states each side of the pinning level, σ 1e-15 cm²). "
   "It depletes its surroundings like a floating Schottky contact and is a strong generation centre "
   "in a depletion region (leakage).", material="Metal (incl.)", pin=True, metal="NiSi",
   shape="sphere", size=(60., 60., 60.))
_p("Cu3Si precipitate (Cu contamination)", "Voids & inclusions", "inclusion",
   "Copper precipitates as metallic Cu3Si particles (Cu barrier on n-Si 0.58 eV): copper's main "
   "electrical form in silicon (Seibt et al., Phys. Status Solidi A 166, 171, 1998). Scatter many with a "
   "density to model a contaminated wafer.",
   material="Metal (incl.)", pin=True, metal="Cu", shape="sphere", size=(30., 30., 30.))
_p("COP (crystal-originated void, Si)", "Voids & inclusions", "inclusion",
   "Octahedral vacancy-cluster void of Czochralski Si (~100–200 nm) with a 2–4 nm oxide on its inner "
   "wall (D_it ~1e12 cm⁻²eV⁻¹); near the surface it thins the gate oxide and degrades its integrity.",
   material="Void", dit=1e12, shape="box", size=(120., 120., 120.))
_p("Micropipe (4H-SiC)", "Voids & inclusions", "inclusion",
   "Hollow-core screw dislocation along the c-axis (y, depth): an open tube 0.5–5 µm across through "
   "the whole epilayer, its wall full of surface states (D_it ~1e13 cm⁻²eV⁻¹). A killer defect for "
   "SiC power devices (leakage, early breakdown). Length 0 = through the whole device.",
   material="Void", dit=1e13, shape="cylinder", size=(500., 0., 0.), axis="y", host="4H-SiC")

# ---- point defects, contamination (Si)
_p("Deep trap (custom)", "Point defects & contamination", "volume",
   "Any deep level: set its type, energy and capture cross sections below. SRH recombination through "
   "the level plus the charge it traps (acceptor: − when filled; donor: + when empty).",
   levels=[L("A", "Ei", 0.0, 1e-15, 1e-15)], density=1e15, shape="sphere", size=(200., 200., 200.))
_p("Fe_i (interstitial iron, Si)", "Point defects & contamination", "uniform",
   "Interstitial iron, the most common lifetime killer in p-type Si: donor Ev + 0.38 eV, "
   "σn 5e-14, σp 7e-17 cm² (Istratov, Hieslmair, Weber, Appl. Phys. A 69, 13, 1999). "
   "Typical contamination 1e10–1e13 cm⁻³.",
   levels=[L("D", "Ev", 0.38, 5e-14, 7e-17)], density=1e12, shape="everywhere", host="Si")
_p("FeB pair (p-Si)", "Point defects & contamination", "uniform",
   "Iron–boron pair (Fe_i paired with B in p-Si at room temperature): acceptor Ec − 0.23 eV, "
   "σn 5e-15, σp 3e-15 cm² (Rein & Glunz, APL 82, 1054, 2003).",
   levels=[L("A", "Ec", 0.23, 5e-15, 3e-15)], density=1e12, shape="everywhere", host="Si")
_p("Au (gold, lifetime control)", "Point defects & contamination", "uniform",
   "Substitutional gold, classic lifetime killer of fast diodes: acceptor Ec − 0.55 eV (σn 8.5e-17, "
   "σp 9e-15 cm²) and donor Ev + 0.35 eV (σn 3.2e-16, σp 3.5e-15 cm²) - representative 300 K values; "
   "measured cross sections scatter (Bullis, Solid-State Electron. 9, 143, 1966; Graff, Metal Impurities "
   "in Silicon-Device Fabrication, 2000). Lifetime control: 1e13–1e15 cm⁻³; it also compensates "
   "lightly doped regions.",
   levels=[L("A", "Ec", 0.55, 8.5e-17, 9e-15), L("D", "Ev", 0.35, 3.2e-16, 3.5e-15)],
   density=1e14, shape="everywhere", host="Si")
_p("Pt (platinum, lifetime control)", "Point defects & contamination", "uniform",
   "Substitutional platinum (power-diode lifetime control, less leakage than Au): acceptor Ec − 0.23 eV, "
   "donor Ev + 0.32 eV. Measured capture cross sections scatter by 2–3 orders of magnitude; "
   "representative values used here: acceptor σn 1e-15/σp 2e-14, donor σn 2e-14/σp 1e-15 cm².",
   levels=[L("A", "Ec", 0.23, 1e-15, 2e-14), L("D", "Ev", 0.32, 2e-14, 1e-15)],
   density=1e14, shape="everywhere", host="Si")
_p("Neutron damage (Si, Perugia model)", "Radiation damage", "fluence",
   "Displacement damage by 1-MeV-equivalent neutrons, the three-level 'Perugia' trap model for "
   "p-type FZ silicon (Petasecca et al., IEEE TNS 53, 2971, 2006): acceptors Ec − 0.42 eV (σn 2e-15, "
   "σp 2e-14, η 1.613 cm⁻¹) and Ec − 0.46 eV (5e-15/5e-14, η 0.9) and a donor Ev + 0.36 eV "
   "(2.5e-14/2.5e-15, η 0.9). Density = fluence Φeq (n_eq/cm²); level densities η·Φ.",
   levels=[L("A", "Ec", 0.42, 2e-15, 2e-14, 1.613), L("A", "Ec", 0.46, 5e-15, 5e-14, 0.9),
           L("D", "Ev", 0.36, 2.5e-14, 2.5e-15, 0.9)], density=1e13, shape="everywhere", host="Si")
_p("Ionizing dose (MOS: oxide charge + D_it)", "Radiation damage", "interface",
   "Total-ionizing-dose damage of a MOS interface: holes trapped in the oxide near the interface "
   "(fixed positive charge, default 5e11 q/cm²) and radiation-generated interface traps "
   "(D_it 5e11 cm⁻²eV⁻¹). Threshold shifts negative, the subthreshold swing grows.",
   dit=5e11, dit_sigma=1e-16, charge=5e11, shape="everywhere")

# ---- interfaces
_p("Interface traps Si/SiO2 (D_it)", "Interfaces", "interface",
   "Interface-state continuum at semiconductor/insulator interfaces inside the shape (and under face "
   "gates): constant D_it across the gap, acceptor-like above midgap, donor-like below, σ 1e-16 cm². "
   "Thermal oxide after anneal ~1e10, as grown or damaged 1e11–1e12 cm⁻²eV⁻¹. Adds C_it = q²D_it to the "
   "subthreshold swing.", dit=1e11, dit_sigma=1e-16, shape="everywhere")
_p("Interface traps 4H-SiC/SiO2", "Interfaces", "interface",
   "SiC/SiO2: D_it ~3e11 cm⁻²eV⁻¹ across the gap plus the near-interface acceptor traps just below Ec "
   "(1e12, 6e11, 3e11 cm⁻² at Ec − 0.1/0.2/0.3 eV) that trap channel electrons and cut the SiC MOSFET's "
   "channel charge.", dit=3e11, dit_sigma=1e-15,
   levels=[L("A", "Ec", 0.10, 1e-15, 1e-15, 1e12), L("A", "Ec", 0.20, 1e-15, 1e-15, 6e11),
           L("A", "Ec", 0.30, 1e-15, 1e-15, 3e11)], shape="everywhere", host="4H-SiC")
_p("Oxide fixed charge (Qf)", "Interfaces", "interface",
   "Fixed charge at semiconductor/insulator interfaces inside the shape (q/cm², + or −): shifts the "
   "flat band and threshold by −Q/C_ox. Good thermal oxide on Si ~5e10, poor 1e11–1e12.",
   charge=1e11, shape="everywhere")
_p("Fixed charge (volume)", "Interfaces", "charge",
   "Fixed space charge (q/cm³, + or −) in the shape - e.g. a charged precipitate cloud.",
   charge=1e16, shape="sphere", size=(200., 200., 200.))

# ---- extended defects
_p("Dislocation (Si)", "Extended defects", "line",
   "Dislocation line with deep states along its core (decorated/contaminated dislocations: "
   "~1e6 states per cm of line, acceptor-like around Ec − 0.45 eV, σ 1e-15 cm²): a recombination-active, "
   "charged line. Cylinder axis = the line; length 0 = through the whole device.",
   levels=[L("A", "Ec", 0.45, 1e-15, 1e-15)], density=1e6, shape="cylinder", size=(10., 0., 0.), axis="y",
   host="Si")
_p("Threading dislocation (GaN)", "Extended defects", "line",
   "Threading dislocation of GaN epitaxy (densities 1e7–1e10 cm⁻²): acceptor-like states along the line "
   "(~1e7 cm⁻¹, Ec − 1.0 eV, σ 1e-16 cm²) - a negatively charged line that scatters and depletes "
   "electrons and conducts leakage. Scatter them with a density per cm².",
   levels=[L("A", "Ec", 1.0, 1e-16, 1e-16)], density=1e7, shape="cylinder", size=(10., 0., 0.), axis="y",
   host="GaN")
_p("Stacking fault (4H-SiC)", "Extended defects", "plane",
   "Basal-plane stacking fault (the bipolar-degradation defect of SiC PiN diodes and MOSFET body diodes): "
   "a recombination-active plane, ~1e12 cm⁻² midgap centres, σ 1e-15 cm². Slab normal = axis; "
   "widths 0 = across the whole device.",
   levels=[L("N", "Ei", 0.0, 1e-15, 1e-15)], density=1e12, shape="slab", size=(2., 0., 0.), axis="y",
   host="4H-SiC")
_p("Grain boundary (poly-Si)", "Extended defects", "plane",
   "Grain boundary of polycrystalline silicon: a sheet of gap states (D_it 5e12 cm⁻²eV⁻¹ across the gap, "
   "amphoteric, σ 1e-15 cm² - ~6e12 cm⁻² states, Seto's 3e12 trapped) that charges up and forms a "
   "potential barrier in each grain boundary, plus recombination.",
   dit=5e12, dit_sigma=1e-15, shape="slab", size=(2., 0., 0.), axis="x", host="Si")

# ---- SiC, GaN, GaAs
_p("Z1/2 (4H-SiC, carbon vacancy)", "SiC · GaN · GaAs", "uniform",
   "The lifetime killer of n-type 4H-SiC: acceptor-like Ec − 0.65 eV, σn 2e-14, σp 1e-14 cm²; "
   "as-grown epilayers 1e12–1e14 cm⁻³ (τ ≈ 1 µs at ~1e13).",
   levels=[L("A", "Ec", 0.65, 2e-14, 1e-14)], density=1e13, shape="everywhere", host="4H-SiC")
_p("EH6/7 (4H-SiC)", "SiC · GaN · GaAs", "uniform",
   "Deep centre of 4H-SiC (carbon vacancy, deeper charge state): donor-like Ec − 1.55 eV, "
   "σn 1e-14, σp 1e-14 cm²; usually tracks Z1/2.",
   levels=[L("D", "Ec", 1.55, 1e-14, 1e-14)], density=1e13, shape="everywhere", host="4H-SiC")
_p("C_N acceptor (GaN, C-doped buffer)", "SiC · GaN · GaAs", "uniform",
   "Carbon on nitrogen site, the deep acceptor of semi-insulating C-doped GaN buffers: Ev + 0.90 eV, "
   "σn 1e-17, σp 1e-15 cm²; buffers 1e17–1e19 cm⁻³ (current collapse, buffer leakage).",
   levels=[L("A", "Ev", 0.90, 1e-17, 1e-15)], density=1e18, shape="everywhere", host="GaN")
_p("Fe_Ga acceptor (GaN, Fe-doped buffer)", "SiC · GaN · GaAs", "uniform",
   "Iron in GaN (semi-insulating buffers of RF HEMTs): acceptor Ec − 0.55 eV, σ 1e-15 cm², 1e17–1e18 cm⁻³.",
   levels=[L("A", "Ec", 0.55, 1e-15, 1e-15)], density=1e17, shape="everywhere", host="GaN")
_p("EL2 (GaAs)", "SiC · GaN · GaAs", "uniform",
   "The native deep donor (As antisite) that makes undoped GaAs semi-insulating: Ec − 0.75 eV, "
   "σn 1e-16, σp 1e-18 cm², ~1.5e16 cm⁻³ in LEC GaAs.",
   levels=[L("D", "Ec", 0.75, 1e-16, 1e-18)], density=1.5e16, shape="everywhere", host="GaAs")

CATEGORIES = []
for _n, _d in LIB.items():
    if _d["cat"] not in CATEGORIES: CATEGORIES.append(_d["cat"])

SHAPES = ["sphere", "box", "cylinder", "slab", "everywhere"]
SIZE_LABELS = {"sphere": ("diameter", "", ""), "box": ("dx", "dy", "dz"),
               "cylinder": ("diameter", "length (0 = all)", ""),
               "slab": ("thickness", "width 1 (0 = all)", "width 2 (0 = all)"),
               "everywhere": ("", "", "")}
DENSITY_UNIT = {"inclusion": "", "volume": "cm⁻³", "uniform": "cm⁻³", "fluence": "n_eq/cm² (fluence)",
                "line": "cm⁻¹ (per length)", "plane": "cm⁻² (per area)", "interface": "", "charge": ""}


def preset_unit(kind):
    return DENSITY_UNIT.get(LIB[kind]["geom"], "") if kind in LIB else ""


# --------------------------------------------------------------- defect
class Defect:
    """One defect.  Geometry in nm; centre (x, y, z); size (a, b, c):
    sphere a = diameter; box a, b, c = dx, dy, dz; cylinder a = diameter,
    b = length along `axis` (0: the whole domain); slab a = thickness along
    the normal `axis`, b, c = widths along the other two axes in xyz order
    (0: whole domain); everywhere: the whole host material.
    density: per the preset's unit; dit: cm^-2 eV^-1 on an inclusion's surface
    or at interfaces; charge: q/cm^2 (interface) or q/cm^3 (charge);
    metal: a metal particle's metal; level: (type, ref, E, sn, sp) override
    of the first level (custom traps); host: material filter ("" = every
    semiconductor)."""
    def __init__(self, kind="Void", shape=None, x=0., y=0., z=0., a=None, b=None, c=None, axis=None,
                 density=None, dit=None, charge=None, metal=None, level=None, host=None, label="", group=""):
        p = LIB.get(kind) or LIB["Void"]
        self.kind = kind if kind in LIB else "Void"
        self.shape = shape or p["shape"]
        self.x, self.y, self.z = float(x), float(y), float(z)
        sz = p["size"]
        self.a = float(sz[0] if a is None else a); self.b = float(sz[1] if b is None else b)
        self.c = float(sz[2] if c is None else c)
        self.axis = axis or p["axis"]
        self.density = float(p["density"] if density is None else density)
        self.dit = float(p["dit"] if dit is None else dit)
        self.charge = float(p["charge"] if charge is None else charge)
        self.metal = p["metal"] if metal is None else metal
        self.level = list(level) if level else None
        self.host = p["host"] if host is None else host
        self.label = label or self.kind
        self.group = group

    # -- persistence
    def to_dict(self):
        return dict(kind=self.kind, shape=self.shape, x=self.x, y=self.y, z=self.z, a=self.a, b=self.b,
                    c=self.c, axis=self.axis, density=self.density, dit=self.dit, charge=self.charge,
                    metal=self.metal, level=self.level, host=self.host, label=self.label, group=self.group)
    @classmethod
    def from_dict(cls, d):
        return cls(**{k: d.get(k) for k in ("kind", "shape", "x", "y", "z", "a", "b", "c", "axis", "density",
                                             "dit", "charge", "metal", "level", "host", "label", "group")
                      if k in d})
    def key(self):
        return repr(sorted(self.to_dict().items()))

    def preset(self):
        return LIB[self.kind]

    def levels(self):
        """The preset's levels, the first replaced by the custom level."""
        lv = [Level(l.type, l.ref, l.E, l.sn, l.sp, l.w) for l in self.preset()["levels"]]
        if self.level and lv:
            t, r, E, sn, sp = self.level
            lv[0] = Level(t, r, E, sn, sp, lv[0].w)
        return lv

    # -- geometry (nm)
    def extent(self, L):
        """Axis-aligned box (lo, hi) the defect occupies in a domain L (nm)."""
        c = np.array([self.x, self.y, self.z], float); Lv = np.asarray(L, float)
        if self.shape == "everywhere":
            return np.zeros(3), Lv.copy()
        if self.shape == "sphere":
            r = 0.5 * self.a; return c - r, c + r
        if self.shape == "box":
            h = 0.5 * np.array([self.a, self.b, self.c]); return c - h, c + h
        ax = AXES.index(self.axis)
        if self.shape == "cylinder":
            r = 0.5 * self.a; lo, hi = c - r, c + r
            if self.b > 0: lo[ax] = c[ax] - 0.5 * self.b; hi[ax] = c[ax] + 0.5 * self.b
            else: lo[ax] = 0.; hi[ax] = Lv[ax]
            return lo, hi
        # slab
        lo, hi = np.zeros(3), Lv.copy()
        lo[ax] = c[ax] - 0.5 * self.a; hi[ax] = c[ax] + 0.5 * self.a
        o = [q for q in range(3) if q != ax]
        for q, w in zip(o, (self.b, self.c)):
            if w > 0: lo[q] = c[q] - 0.5 * w; hi[q] = c[q] + 0.5 * w
        return lo, hi

    def distance(self, P, L):
        """Signed-like distance (nm) of points P (n,3) from the defect's core
        geometry for the line/plane kinds: the distance to the cylinder axis
        segment, or to the slab mid-plane (inf outside the in-plane extent)."""
        lo, hi = self.extent(L); ax = AXES.index(self.axis)
        c = np.array([self.x, self.y, self.z], float)
        if self.shape == "cylinder":
            o = [q for q in range(3) if q != ax]
            d = np.hypot(P[:, o[0]] - c[o[0]], P[:, o[1]] - c[o[1]])
            out = (P[:, ax] < lo[ax] - 1e-9) | (P[:, ax] > hi[ax] + 1e-9)
            return np.where(out, np.inf, d)
        if self.shape == "slab":
            d = np.abs(P[:, ax] - c[ax])
            o = [q for q in range(3) if q != ax]
            out = np.zeros(len(P), bool)
            for q in o: out |= (P[:, q] < lo[q] - 1e-9) | (P[:, q] > hi[q] + 1e-9)
            return np.where(out, np.inf, d)
        return np.full(len(P), np.inf)

    def contains(self, P, L, grow=0.):
        """Points P (n,3, nm) inside the shape grown by `grow` nm."""
        c = np.array([self.x, self.y, self.z], float)
        if self.shape == "everywhere":
            return np.ones(len(P), bool)
        if self.shape == "sphere":
            return np.sum((P - c) ** 2, axis=1) <= (0.5 * self.a + grow) ** 2 * (1 + 1e-12)
        if self.shape == "box":
            # half-open like the regions: a node ON the low face is inside, ON the
            # high face outside (with pinned faces the box is then exact)
            h = 0.5 * np.array([self.a, self.b, self.c]) + grow
            e = 1e-9 * max(float(np.max(h)), 1.)
            return np.all((P >= c - h - e) & (P < c + h - e), axis=1)
        if self.shape == "cylinder":
            return self.distance(P, L) <= 0.5 * self.a + grow
        return self.distance(P, L) <= 0.5 * self.a + grow

    def measure(self, L):
        """Volume (nm^3), length (nm) or area (nm^2) of the core geometry."""
        lo, hi = self.extent(L)
        if self.shape == "sphere": return math.pi / 6. * self.a ** 3
        if self.shape == "box": return self.a * self.b * self.c
        if self.shape == "everywhere": return float(np.prod(hi - lo))
        ax = AXES.index(self.axis)
        if self.shape == "cylinder":
            return math.pi * 0.25 * self.a ** 2 * (hi[ax] - lo[ax])
        o = [q for q in range(3) if q != ax]
        return (hi[o[0]] - lo[o[0]]) * (hi[o[1]] - lo[o[1]]) * self.a

    def line_length(self, L):
        lo, hi = self.extent(L); ax = AXES.index(self.axis); return hi[ax] - lo[ax]

    def plane_area(self, L):
        lo, hi = self.extent(L); ax = AXES.index(self.axis)
        o = [q for q in range(3) if q != ax]
        return (hi[o[0]] - lo[o[0]]) * (hi[o[1]] - lo[o[1]])

    def describe(self):
        g = self.preset()["geom"]
        s = {"sphere": f"⌀{self.a:g}", "box": f"{self.a:g}×{self.b:g}×{self.c:g}",
             "cylinder": f"⌀{self.a:g} ∥{self.axis}" + (f" L{self.b:g}" if self.b > 0 else ""),
             "slab": f"t{self.a:g} ⊥{self.axis}", "everywhere": "everywhere"}[self.shape]
        return s

    def strength_text(self):
        g = self.preset()["geom"]
        if g in ("volume", "uniform"): return f"{self.density:.2g} cm⁻³"
        if g == "fluence": return f"Φ {self.density:.2g} cm⁻²"
        if g == "line": return f"{self.density:.2g} cm⁻¹"
        if g == "plane": return (f"{self.density:.2g} cm⁻²" if self.preset()["levels"] else f"D_it {self.dit:.2g}")
        if g == "inclusion": return (f"{self.metal}" if self.preset()["pin"] else f"D_it {self.dit:.2g}")
        if g == "interface":
            t = []
            if self.dit: t.append(f"D_it {self.dit:.2g}")
            if self.charge: t.append(f"Q {self.charge:+.2g}")
            return ", ".join(t) or "-"
        if g == "charge": return f"{self.charge:+.2g} q/cm³"
        return ""


def scatter(proto, n, box, seed=1, spread=0.0, group=None):
    """n copies of the defect `proto` at uniformly random centres inside the
    box (lo, hi) (nm), sizes varied by +-spread (fraction), reproducible
    with `seed`."""
    rng = np.random.default_rng(int(seed))
    lo, hi = np.asarray(box[0], float), np.asarray(box[1], float)
    out = []
    grp = group or f"{proto.kind} ×{n} (seed {seed})"
    for k in range(int(n)):
        d = Defect.from_dict(proto.to_dict())
        P = lo + rng.random(3) * (hi - lo)
        d.x, d.y, d.z = (float(v) for v in P)
        if spread > 0:
            f = 1. + spread * (2. * rng.random() - 1.)
            d.a *= f
            if d.shape == "box": d.b *= f; d.c *= f
        d.label = f"{proto.label} #{k+1}"; d.group = grp
        out.append(d)
    return out


def count_for_density(proto, dens, L, box=None):
    """Number of defects for a density: per cm^2 of the plane normal to the
    defect's axis for lines (e.g. a threading-dislocation density), per cm^3
    of the box otherwise."""
    lo, hi = (np.zeros(3), np.asarray(L, float)) if box is None else (np.asarray(box[0], float), np.asarray(box[1], float))
    ext = (hi - lo) * 1e-7                          # cm
    if proto.shape == "cylinder" and proto.b <= 0:
        ax = AXES.index(proto.axis); o = [q for q in range(3) if q != ax]
        return dens * ext[o[0]] * ext[o[1]]
    return dens * float(np.prod(ext))


# --------------------------------------------------------------- meshes
class MeshView:
    """Node positions (nm), volumes (m^3), per-node material names and the
    edges (with face areas, m^2) of a tensor or prism mesh."""
    def __init__(self, mat_names, tensor=None, graph=None):
        self.mat = np.asarray(mat_names, dtype=object)
        if tensor is not None:
            xs, ys, zs = (np.asarray(a, float) for a in tensor)
            self.kind = "tensor"; self.c = (xs, ys, zs)
            self.n = (len(xs), len(ys), len(zs)); self.N = len(xs) * len(ys) * len(zs)
            self.fv = tuple(self._fv(a) for a in self.c)
        else:
            N, rp, cj, ar, vol, pos = graph
            self.kind = "graph"; self.N = int(N)
            self.rp = np.asarray(rp, np.int64); self.cj = np.asarray(cj, np.int64)
            self.ar = np.asarray(ar, float); self.vol = np.asarray(vol, float)
            self.pos = np.asarray(pos, float).reshape(-1, 3) * 1e9
            self._grid = None

    @staticmethod
    def _fv(c):
        if len(c) < 2: return np.ones(len(c))
        h = np.diff(c); f = np.empty(len(c)); f[0] = 0.5 * h[0]; f[-1] = 0.5 * h[-1]
        f[1:-1] = 0.5 * (h[1:] + h[:-1]); return f

    # -- queries
    def box_nodes(self, lo, hi):
        """Indices of the nodes inside [lo, hi] (nm)."""
        lo = np.asarray(lo, float) * 1e-9; hi = np.asarray(hi, float) * 1e-9
        if self.kind == "tensor":
            r = []
            for q in range(3):
                c = self.c[q]
                a = int(np.searchsorted(c, lo[q] - 1e-15, "left")); b = int(np.searchsorted(c, hi[q] + 1e-15, "right"))
                r.append(np.arange(a, b))
            if not all(len(x) for x in r): return np.zeros(0, np.int64)
            I, J, K = np.meshgrid(r[0], r[1], r[2], indexing="ij")
            return (I + self.n[0] * (J + self.n[1] * K)).ravel().astype(np.int64)
        lo9, hi9 = lo * 1e9 - 1e-6, hi * 1e9 + 1e-6
        if self._grid is None:                      # coarse buckets, built once
            P = self.pos; p0 = P.min(0); p1 = P.max(0); nb = 24
            w = np.maximum((p1 - p0) / nb, 1e-9)
            b = np.clip(((P - p0) / w).astype(np.int64), 0, nb - 1)
            key = b[:, 0] + nb * (b[:, 1] + nb * b[:, 2])
            order = np.argsort(key, kind="stable"); ks = key[order]
            starts = np.searchsorted(ks, np.arange(nb ** 3 + 1))
            self._grid = (p0, w, nb, order, starts)
        p0, w, nb, order, starts = self._grid
        b0 = np.clip(((lo9 - p0) / w).astype(np.int64), 0, nb - 1); b1 = np.clip(((hi9 - p0) / w).astype(np.int64), 0, nb - 1)
        cand = []
        for kz in range(b0[2], b1[2] + 1):
            for ky in range(b0[1], b1[1] + 1):
                a = b0[0] + nb * (ky + nb * kz); e = b1[0] + nb * (ky + nb * kz)
                cand.append(order[starts[a]:starts[e + 1]])
        if not cand: return np.zeros(0, np.int64)
        c = np.concatenate(cand)
        m = np.all((self.pos[c] >= lo9) & (self.pos[c] <= hi9), axis=1)
        return np.sort(c[m]).astype(np.int64)

    def positions(self, idx):
        if self.kind == "tensor":
            i = idx % self.n[0]; j = (idx // self.n[0]) % self.n[1]; l = idx // (self.n[0] * self.n[1])
            return np.stack([self.c[0][i], self.c[1][j], self.c[2][l]], 1) * 1e9
        return self.pos[idx]

    def volumes(self, idx):
        if self.kind == "tensor":
            i = idx % self.n[0]; j = (idx // self.n[0]) % self.n[1]; l = idx // (self.n[0] * self.n[1])
            return self.fv[0][i] * self.fv[1][j] * self.fv[2][l]
        return self.vol[idx]

    def spacing(self, idx):
        """Local node spacing (nm): the cube root of the control volume."""
        return np.cbrt(np.maximum(self.volumes(idx), 1e-60)) * 1e9

    def edges(self, idx):
        """(k, kb, area m^2) of every edge from the nodes idx."""
        idx = np.asarray(idx, np.int64)
        if self.kind == "graph":
            if not len(idx): return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0)
            st = self.rp[idx]; cnt = self.rp[idx + 1] - st
            k = np.repeat(idx, cnt)
            first = np.concatenate(([0], np.cumsum(cnt)[:-1]))
            off = np.repeat(st - first, cnt) + np.arange(int(cnt.sum()))
            return k, self.cj[off], self.ar[off]
        nx, ny, nz = self.n
        i = idx % nx; j = (idx // nx) % ny; l = idx // (nx * ny)
        K, KB, A = [], [], []
        for ax, (st, nn) in enumerate(((1, nx), (nx, ny), (nx * ny, nz))):
            pos = (i, j, l)[ax]
            o = [q for q in range(3) if q != ax]
            area = self.fv[o[0]][(i, j, l)[o[0]]] * self.fv[o[1]][(i, j, l)[o[1]]]
            for sg in (-1, 1):
                ok = (pos + sg >= 0) & (pos + sg < nn)
                K.append(idx[ok]); KB.append(idx[ok] + sg * st); A.append(area[ok])
        return np.concatenate(K), np.concatenate(KB), np.concatenate(A)


# --------------------------------------------------------------- mapping
class DefectData:
    """What the defects do to a mesh: material overrides {node: name}, deep
    levels [(level key (type, ref, E, sn_m2, sp_m2), nodes, densities m^-3)],
    fixed charge (nodes, q/m^3), and notes."""
    def __init__(self):
        self.mat = {}; self.traps = []; self.qn = []; self.qv = []; self.notes = []
        self._agg = {}                         # (who, message) -> [count, values], per group
    def note(self, d, msg, val=None):
        """A note about defect d; the members of a scattered group are counted
        together ('{}' in msg: the range of the values val)."""
        who = f"group '{d.group}'" if d.group else f"'{d.label}'"
        e = self._agg.setdefault((who, msg), [0, []]); e[0] += 1
        if val is not None: e[1].append(val)
    def finish(self):
        for (who, msg), (n, vals) in self._agg.items():
            if vals:
                lo, hi = min(vals), max(vals)
                msg = msg.replace("{}", f"{lo:.0f}" if hi - lo < 0.5 else f"{lo:.0f}-{hi:.0f}")
            cnt = f" ({n} defect{'s' if n > 1 else ''})" if who.startswith("group") else ""
            if n > 1: msg = msg.replace(" its ", " their ")
            self.notes.append(f"{who}{cnt}: {msg}")
        self._agg = {}
    def add_level(self, lv, nodes, dens):
        nodes = np.asarray(nodes, np.int64); dens = np.asarray(dens, float)
        ok = dens > 0
        if not ok.any(): return
        key = (lv.type, lv.ref, float(lv.E), lv.sn * 1e-4, lv.sp * 1e-4)
        self.traps.append((key, nodes[ok], dens[ok]))
    def add_charge(self, nodes, rho):
        self.qn.append(np.asarray(nodes, np.int64)); self.qv.append(np.asarray(rho, float))
    def charge(self):
        if not self.qn: return np.zeros(0, np.int64), np.zeros(0)
        return np.concatenate(self.qn), np.concatenate(self.qv)
    def n_entries(self):
        return sum(len(t[1]) for t in self.traps)


def dit_levels(dit, sigma, Eg, nb=24):
    """A constant D_it (cm^-2 eV^-1) across a gap Eg as nb amphoteric levels
    (acceptor above midgap, donor below): [(Level, sheet density cm^-2)]."""
    dE = Eg / nb; out = []
    for k in range(nb):
        E = (k + 0.5) * dE                       # above Ev
        t = "A" if E > 0.5 * Eg else "D"
        out.append((Level(t, "Ev", E, sigma, sigma), dit * dE))
    return out


def _fib_sphere(n):
    """n nearly uniform unit vectors (Fibonacci lattice)."""
    i = np.arange(n) + 0.5
    ph = np.arccos(1. - 2. * i / n); th = math.pi * (1. + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)], 1)


def inclusion_geometry(d, L):
    """(volume nm^3, surface area nm^2) of an inclusion's shape inside the
    domain (0, L): exact for shapes inside it and for boxes, sampled for the
    clipped parts of spheres and cylinders.  Surface on a domain face does
    not count (it is not an interface with the semiconductor)."""
    Lv = np.asarray(L, float); c = np.array([d.x, d.y, d.z], float)
    lo, hi = d.extent(L); tol = 1e-9 * max(float(Lv.max()), 1.)
    inn = lambda P: np.all((P > tol) & (P < Lv - tol), axis=1)
    if d.shape == "box":
        lc = np.maximum(lo, 0.); hc = np.minimum(hi, Lv); e = np.maximum(hc - lc, 0.)
        A = 0.
        for q in range(3):
            o = [k for k in range(3) if k != q]
            for v in (lo[q], hi[q]):
                if tol < v < Lv[q] - tol: A += e[o[0]] * e[o[1]]
        return float(np.prod(e)), float(A)
    def vol_sampled():
        lc = np.maximum(lo, 0.); hc = np.minimum(hi, Lv)
        if np.any(hc <= lc): return 0.
        n = 40; g = [lc[q] + (np.arange(n) + 0.5) / n * (hc[q] - lc[q]) for q in range(3)]
        X, Y, Z = np.meshgrid(*g, indexing="ij")
        P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
        return float(d.contains(P, L).mean() * np.prod(hc - lc))
    r = 0.5 * d.a
    if d.shape == "sphere":
        V0, A0 = math.pi / 6. * d.a ** 3, math.pi * d.a ** 2
        if np.all(lo > tol) and np.all(hi < Lv - tol): return V0, A0
        return vol_sampled(), A0 * float(inn(c + r * _fib_sphere(4000)).mean())
    # cylinder along its axis q: lateral surface + the caps of a finite one
    q = AXES.index(d.axis); o = [k for k in range(3) if k != q]
    l0, l1 = max(lo[q], 0.), min(hi[q], Lv[q]); ln = max(l1 - l0, 0.)
    side_in = all(c[k] - r > tol and c[k] + r < Lv[k] - tol for k in o)
    caps = [v for v in (lo[q], hi[q]) if tol < v < Lv[q] - tol]
    if side_in:
        return math.pi * r * r * ln, 2. * math.pi * r * ln + math.pi * r * r * len(caps)
    th = np.linspace(0., 2. * math.pi, 361)[:-1]; sv = l0 + (np.arange(60) + 0.5) / 60. * ln
    T, S = np.meshgrid(th, sv, indexing="ij")
    P = np.zeros((T.size, 3)); P[:, q] = S.ravel()
    P[:, o[0]] = c[o[0]] + r * np.cos(T.ravel()); P[:, o[1]] = c[o[1]] + r * np.sin(T.ravel())
    A = 2. * math.pi * r * ln * float(inn(P).mean())
    k = np.arange(800) + 0.5; rr = r * np.sqrt(k / 800.); aa = math.pi * (3. - 5 ** 0.5) * k
    for v in caps:
        D = np.zeros((800, 3)); D[:, q] = v
        D[:, o[0]] = c[o[0]] + rr * np.cos(aa); D[:, o[1]] = c[o[1]] + rr * np.sin(aa)
        A += math.pi * r * r * float(inn(D).mean())
    return vol_sampled(), A


def _nearest_host(view, d, ok_mask, max_pad=1e4):
    """The host node nearest to a defect's centre (array of one index), or empty."""
    c = np.array([d.x, d.y, d.z], float); h = 5.
    while h <= max_pad:
        idx = view.box_nodes(c - h, c + h)
        if len(idx):
            idx = idx[ok_mask(idx)]
            if len(idx):
                dist = np.sum((view.positions(idx) - c) ** 2, axis=1)
                return idx[[int(np.argmin(dist))]]
        h *= 2.
    return np.zeros(0, np.int64)


def _spread(dd, view, idx, total, P=None):
    """Spread `total` (a number of states) uniformly over the nodes idx:
    per-node density (m^-3) so that sum(N_k V_k) = total."""
    V = view.volumes(idx); s = V.sum()
    return np.full(len(idx), total / s) if s > 0 else np.zeros(len(idx))


def apply(defects, view, L, semis, insulators, eg_of, barrier_of=None, gate_faces=()):
    """Map defects onto a mesh.  view: MeshView; L: domain (nm); semis:
    material names that are semiconductors; insulators: names of insulators
    (voids, oxides, metal inclusions); eg_of(name) -> band gap (eV);
    barrier_of(metal, semi) -> phi_Bn (eV); gate_faces: [(nodes, areas m^2)]
    of the semiconductor nodes under face (oxide) gates - interfaces too.
    Returns DefectData."""
    out = DefectData()
    mat = view.mat.copy()                              # running material map
    issemi = np.isin(mat, np.array(sorted(semis), dtype=object))
    for d in defects:
        p = d.preset(); g = p["geom"]
        lo, hi = d.extent(L)
        def host_ok(names, _h=d.host):          # names: mat[idx] - a host-material filter
            return (names == _h) if _h else np.ones(len(names), bool)
        if g == "inclusion":
            idx = view.box_nodes(lo, hi)
            P = view.positions(idx)
            ins = idx[d.contains(P, L) & host_ok(mat[idx]) & issemi[idx]]
            Vs, As = inclusion_geometry(d, L)                  # nm^3, nm^2 inside the domain
            if Vs <= 0.: out.note(d, "outside the device"); continue
            Vin = float(view.volumes(ins).sum()) * 1e27 if len(ins) else 0.
            def states(uk, area_m2, Vk, name):
                """Surface states (per level: m^-3 at the nodes uk) for a surface area."""
                if p["pin"]:
                    pb = barrier_of(d.metal, name) if (barrier_of and d.metal) else 0.5 * eg_of(name)
                    Ns = 5e13 * 1e4                                    # m^-2 per level
                    out.add_level(Level("D", "Ec", pb + 0.05, 1e-15, 1e-15), uk, Ns * area_m2 / Vk)
                    out.add_level(Level("A", "Ec", max(pb - 0.05, 0.01), 1e-15, 1e-15), uk, Ns * area_m2 / Vk)
                elif d.dit > 0:
                    for lv, ns in dit_levels(d.dit, p["dit_sigma"], eg_of(name)):
                        out.add_level(lv, uk, ns * 1e4 * area_m2 / Vk)
            if not len(ins) or Vin > 3. * Vs:
                # smaller than the mesh there: replacing whole control volumes would
                # make it several times its size - it acts through its surface
                # states only (exact number), on the nodes inside it or the nearest one
                nod = ins if len(ins) else _nearest_host(view, d, lambda q: issemi[q] & host_ok(mat[q]))
                if not len(nod): out.note(d, "no host material near it"); continue
                Vall = float(view.volumes(nod).sum())
                for name in set(mat[nod]):
                    m = mat[nod] == name; Vk = view.volumes(nod[m])
                    states(nod[m], As * 1e-18 * Vk / Vall, Vk, name)     # area shared by volume
                out.note(d, "smaller than the mesh cells there - modelled by its surface states only "
                            "(refine the mesh to resolve its volume)")
                continue
            if not 0.5 <= Vin / Vs <= 2.:
                out.note(d, "meshed with {} % of its volume - refine the mesh near it", 100. * Vin / Vs)
            for k in ins: out.mat[int(k)] = p["material"]
            mat[ins] = p["material"]; issemi[ins] = False
            # surface states on the semiconductor side of its surface: the staircase
            # of mesh faces carries the TRUE surface area (a meshed sphere's faces
            # add up to ~1.5x its area), so the number of states is exact
            k, kb, A = view.edges(ins)
            outside = ~np.isin(kb, ins)
            Astair = float(A[outside].sum())
            sel = outside & issemi[kb]
            sk = kb[sel]; sa = A[sel] * (As * 1e-18 / Astair if Astair > 0 else 1.)
            if not len(sk): continue
            uk, inv = np.unique(sk, return_inverse=True)
            area = np.bincount(inv, weights=sa)
            Vk = view.volumes(uk)
            for name in set(mat[uk]):
                m = mat[uk] == name
                states(uk[m], area[m], Vk[m], name)
        elif g in ("volume", "uniform", "fluence", "charge"):
            if d.shape == "everywhere":
                idx = np.arange(view.N, dtype=np.int64)
                if g == "charge": nod = idx[issemi[idx]]
                else: nod = idx[issemi[idx] & host_ok(mat[idx])]
                dens_m3 = (d.charge if g == "charge" else d.density) * 1e6
                if not len(nod): out.note(d, "no host material in the device"); continue
                if g == "charge": out.add_charge(nod, np.full(len(nod), dens_m3)); continue
                for lv in d.levels():
                    w = lv.w if g == "fluence" else lv.w
                    out.add_level(lv, nod, np.full(len(nod), w * dens_m3))
                continue
            # a finite shape: keep the number of states (or charges) exact
            grow = 0.
            idx = view.box_nodes(lo, hi); P = view.positions(idx)
            ins = idx[d.contains(P, L) & (issemi[idx] if g == "charge" else host_ok(mat[idx]) & issemi[idx])]
            Vs = d.measure(L) * 1e-27
            Vin = view.volumes(ins).sum() if len(ins) else 0.
            if not (0.5 * Vs <= Vin <= 2. * Vs):          # poorly resolved: spread over the nearest cells
                c = np.array([d.x, d.y, d.z]); h = 0.
                for _ in range(6):
                    h = max(h * 2., 1.)
                    idx2 = view.box_nodes(lo - h, hi + h)
                    if not len(idx2): continue
                    hl = 0.75 * float(np.max(view.spacing(idx2)))
                    idx3 = view.box_nodes(lo - hl, hi + hl); P3 = view.positions(idx3)
                    m3 = d.contains(P3, L, grow=hl) & (issemi[idx3] if g == "charge" else host_ok(mat[idx3]) & issemi[idx3])
                    if m3.any(): ins = idx3[m3]; break
                if not len(ins):
                    out.note(d, "no host node near it"); continue
            amount = (d.charge if g == "charge" else d.density) * 1e6 * Vs     # charges or states
            dens = _spread(out, view, ins, amount)
            if g == "charge": out.add_charge(ins, dens); continue
            for lv in d.levels():
                out.add_level(lv, ins, lv.w * dens)
        elif g in ("line", "plane"):
            # the nodes on the line/plane: those within its radius/half-thickness,
            # or (coarse mesh) the nearest row of host nodes - the number of
            # states is kept exact either way
            nod = np.zeros(0, np.int64)
            for pad in (0., 5., 25., 100., 500., 2000., 1e4):
                idx = view.box_nodes(lo - pad, hi + pad)
                if not len(idx): continue
                P = view.positions(idx); dist = d.distance(P, L)
                hostm = host_ok(mat[idx]) & issemi[idx] & np.isfinite(dist)
                if not hostm.any(): continue
                reach = max(0.5 * d.a, float(dist[hostm].min()) * (1 + 1e-6) + 1e-6)
                if reach > 0.5 * d.a + pad + 1e-9: continue   # a nearer row may lie outside this box
                nod = idx[hostm & (dist <= reach)]
                break
            if not len(nod): out.note(d, "no host node on it"); continue
            if g == "line":
                total = d.density * d.line_length(L) * 1e-7          # states (cm^-1 x cm)
                dens = _spread(out, view, nod, total)
                for lv in d.levels(): out.add_level(lv, nod, lv.w * dens)
            else:
                area_cm2 = d.plane_area(L) * 1e-14
                if p["levels"]:
                    dens = _spread(out, view, nod, d.density * area_cm2)
                    for lv in d.levels(): out.add_level(lv, nod, lv.w * dens)
                if d.dit > 0:
                    for name in set(mat[nod]):
                        m = mat[nod] == name
                        for lv, ns in dit_levels(d.dit, p["dit_sigma"], eg_of(name)):
                            out.add_level(lv, nod[m], _spread(out, view, nod[m], ns * area_cm2 * float(m.sum()) / len(nod)))
        elif g == "interface":
            if d.shape == "everywhere": cand = np.arange(view.N, dtype=np.int64)
            else:
                cand = view.box_nodes(lo, hi); P = view.positions(cand); cand = cand[d.contains(P, L)]
            sem = cand[issemi[cand] & host_ok(mat[cand])]
            k, kb, A = view.edges(sem)
            isin = np.isin(mat[kb], np.array(sorted(insulators), dtype=object)) if len(kb) else np.zeros(0, bool)
            inside = np.isin(kb, cand) if d.shape != "everywhere" else np.ones(len(kb), bool)
            sel = isin & inside
            sk, sa, ok_ = k[sel], A[sel], kb[sel]
            pairs = [(sk, sa, ok_)]
            for gn, ga in gate_faces:              # oxide gates on a domain face
                gn = np.asarray(gn, np.int64); ga = np.asarray(ga, float)
                m = issemi[gn] & host_ok(mat[gn]) & (np.isin(gn, cand) if d.shape != "everywhere" else True)
                pairs.append((gn[m], ga[m], None))
            for sk, sa, ok_ in pairs:
                if not len(sk): continue
                uk, inv = np.unique(sk, return_inverse=True)
                area = np.bincount(inv, weights=sa); Vk = view.volumes(uk)
                for name in set(mat[uk]):
                    m = mat[uk] == name
                    if d.dit > 0:
                        for lv, ns in dit_levels(d.dit, p["dit_sigma"], eg_of(name)):
                            out.add_level(lv, uk[m], ns * 1e4 * area[m] / Vk[m])
                    for lv in d.levels():                   # discrete interface levels (cm^-2 via w)
                        out.add_level(lv, uk[m], lv.w * 1e4 * area[m] / Vk[m])
                if d.charge:
                    if ok_ is not None and len(ok_):        # meshed oxide: on its interface nodes
                        uo, inv2 = np.unique(ok_, return_inverse=True)
                        ao = np.bincount(inv2, weights=sa)
                        out.add_charge(uo, d.charge * 1e4 * ao / view.volumes(uo))
                    else:                                    # face gate: at the semiconductor surface
                        out.add_charge(uk, d.charge * 1e4 * area / Vk)
            if not any(len(s[0]) for s in pairs):
                out.note(d, "no semiconductor/insulator interface inside it")
    out.finish()
    return out


REFINE_MAX = 40          # defects that get their own mesh planes (random scatter: the first ones)


def pin_planes(defects, L):
    """Faces (nm) of box-shaped inclusions (voids, particles), per axis, for
    the mesher's node pairs that make the box exact on a tensor mesh."""
    out = [set(), set(), set()]; n = 0
    for d in defects:
        if d.preset()["geom"] != "inclusion" or d.shape != "box": continue
        n += 1
        if n > REFINE_MAX: break
        lo, hi = d.extent(L)
        for q in range(3):
            for v in (lo[q], hi[q]):
                if 0. < v < L[q]: out[q].add(round(float(v), 6))
    return [sorted(s) for s in out]


def refine_planes(defects, L):
    """Dense planes (nm) per axis that resolve material-changing defects and
    lines/planes on a graded tensor mesh: [xs, ys, zs]."""
    out = [set(), set(), set()]; n = 0
    for d in defects:
        g = d.preset()["geom"]
        if g not in ("inclusion", "line", "plane") or d.shape == "everywhere": continue
        n += 1
        if n > REFINE_MAX: break
        lo, hi = d.extent(L)
        for q in range(3):
            for v in (lo[q], hi[q], 0.5 * (lo[q] + hi[q])):
                if 0. < v < L[q]: out[q].add(round(float(v), 6))
    return [sorted(s) for s in out]
