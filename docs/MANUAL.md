# SemiSim 3D 0.2 — user manual

The long form of the program's Help (F1). Short version: [README](../README.md);
every reference: [SOURCES](SOURCES.md).

**Status.** 0.2 is an early public version and it is still being tested. The
physics is checked against textbook theory (see Validation), not yet against
measurements of real devices; the effects under Limitations are not modelled.
Check important results independently. Free to use, copy, change and share
under the MIT License ([LICENSE](../LICENSE)); no warranty.

SemiSim 3D solves Poisson's equation and the electron and hole continuity
equations (drift-diffusion, box method, Scharfetter-Gummel fluxes) in 3-D. A
device is built from blocks of material and doping - boxes, or boxes cut to a
round or sloped outline - on one of two meshes: a graded rectangular grid
(refined at junctions, gates and interfaces, and on request adapted to the
solution), or triangular prisms (`semimesh.py`) that follow round and sloped
interfaces. The core (`semiconductor_3d.c`, OpenMP) returns bands, quasi-Fermi
levels, fields, currents, mobilities, recombination and trapped charge; the
interface (`main.py`, Python/Tk) builds the device, runs bias points and sweeps
and draws the results. The 55 templates range from textbook diodes to FinFETs,
a CMOS inverter, SiC/GaN power devices, comb-emitter power BJTs and parts
modelled on real components (device class only - not their real dies).

## Install

Windows (WSL2 with Ubuntu) and Linux, in a terminal:

```bash
git clone https://github.com/kotmariusz/Semiconductor_devices_simulator.git
cd Semiconductor_devices_simulator
bash install.sh
```

macOS: first `xcode-select --install` (the compiler) and Python from
[python.org](https://www.python.org/downloads/macos/) or Homebrew, then the
same three commands in Terminal.

`install.sh` installs what is missing (Python 3.8+ with Tk 8.6+, NumPy,
Matplotlib, a C compiler - on Linux from the distribution, with sudo; on macOS
in a virtual environment `.venv`, plus Homebrew's `libomp` for multi-core
solving), builds the solver core (`semibuild.py`), runs a short self-test and
adds the command `semisim3d` - on Linux also a menu entry (on WSL2 in the
Windows Start menu), on macOS `SemiSim3D.command` to double-click. Options:
`PYTHON=/path/to/python3`, `SEMISIM_VENV=1` (Linux: pip instead of system
packages), `NO_OPENMP=1` (macOS: no libomp).

Start: `semisim3d`, or `python3 main.py` in the folder. Update: `git pull` -
the core rebuilds itself on the next start.

By hand: `python3 -m pip install -r requirements.txt`, `python3 semibuild.py`
(`--portable` for any CPU of this kind, `--check` to see the compiler and
OpenMP in use), `python3 main.py`.

Notes: native Windows is not supported (use WSL2). macOS's own
`/usr/bin/python3` (Tk 8.5) does not work; Tk 9 needs Matplotlib 3.10+. Without
OpenMP the core runs on one core. The solver uses all cores (on Apple silicon
the performance cores); Solver page → CPU threads. A core built for one CPU
should not be copied to another computer (`main.py` rebuilds it when needed).

## Using the GUI

The window adapts to the screen: 1280×720 and larger work, and the layout is
remembered between sessions (in `~/.semisim3d.json`). F11 (View → Fill the
screen / restore) makes the window as large as the screen and back; View →
Reset window size brings back the default.

**WSL2 (WSLg) and maximised windows.** WSLg sends the mouse clicks of a
maximised Linux window to the wrong place: after the maximise button the
buttons stop responding, while a drag on a plot still turns the 3-D view.
A normal window as large as the screen does not have the problem, so on WSL
SemiSim turns a maximise (the button, a double click on the title bar,
Win+↑) into F11: the window is un-maximised and given the same rectangle.
View → `WSL: turn maximise into fill screen` switches this off.

```
┌ File  Devices  Run  View  Help ──────────────────────────────────────────────┐
│ ☰  Device [CMOS Inverter ▾]  ▶ Run  ↔ Sweep  ■ Stop  [████░░] point 7/26…   ? │
├──────────── editor (drawer) ───────┬── plots ─────────────────────────────────┤
│ Device│Regions│Contacts│Defects│   │ Structure│3D Slicer│Bands│Current│I-V│… │
│ Solver│Sweep                       │                                          │
│  (every page scrolls both ways;    │   (matplotlib toolbar under every plot)  │
│   wheel / Shift+wheel)             │                                          │
└────────────────────────────────────┴──────────────────────────────────────────┘
```

**Toolbar**

* `☰` (Ctrl+E; ⌘E on a Mac) hides the editor so the plots get the whole window. Drag the
  divider to resize the editor.
* `Device ▾` opens the template library by category. Help → Template library
  shows the full description of each device.
* `▶ Run` (F5) solves at the contact voltages. `↔ Sweep` (F6) runs the sweep.
  `■ Stop` (Esc) stops within one Gummel iteration and keeps the last
  converged bias point.
* The progress bar tracks the bias ramp or the sweep points. Beside it are
  the continuation step, the Gummel iteration and the residual.

**Editor pages**

* `Device`: the template description, domain, base mesh, the mesh type
  (`Rectangular` or `Triangular (prisms)`, with the prism axis), adaptive
  mesh refinement and the final mesh statistics.
* `Regions` and `Contacts`: scrollable tables. Select a row to edit it. A
  contact can be made of a metal (`Metal`; see below).
* `Defects`: the defect table, the defect form (library preset, shape, size,
  density, levels) and random scatter (see Defects).
* `Solver`: temperature, physical models (band-gap narrowing, τ(N), velocity
  saturation, quantum confinement, Lombardi surface mobility, gate
  tunnelling, metal-contact thermionic emission and barrier tunnelling,
  impact ionisation), the background SRH lifetime, iteration limits,
  tolerances and CPU threads. Each template sets the model switches it was
  designed with.
* `Sweep`: the swept terminal, range and points, a floating terminal and the
  plotted current.

**Plot tabs**

* Structure (3-D view plus an XY/YZ/XZ cross-section with the real mesh -
  grid lines, or the triangles of a prism mesh; outlines drawn as they are;
  defects filled where the plane cuts them and dotted otherwise; hover to
  read the region and defect under the cursor), 3-D Slicer, Bands, Current,
  I-V, and Field & Mobility.
* `Results`: every terminal's applied and solved voltage, current, current
  density and resolution, the parameters extracted from the last sweep, and a
  run log.
* `Current`: the number boxes next to `2-D arrows` and `3-D arrows` set how
  many current arrows are drawn - in the in-plane map the arrows across it
  (4-120, the rows follow the map's shape), in the 3-D field the arrows along
  x (2-40; depth y gets 3/4 and z 5/8 of them, never more than the mesh
  lines along an axis, at most 4000 arrows). The current is interpolated
  between the mesh nodes, so the count is not limited by the mesh; denser
  arrows are drawn shorter and thinner. Arrow length is log-scaled, so weak
  and strong currents both show.
* `Current` → `Map` puts the device structure under the arrows of the
  in-plane map: `|J|` (the current density, as before), `structure` (the
  region colours) or `|J| + structure` (both, the structure see-through on
  top). `Structure in 3-D` draws the regions and contacts around the 3-D
  arrows. The `Opacity` slider next to them sets how see-through that
  structure is; parts hidden in `Layers…` stay hidden there, and a part given
  its own opacity in `Layers…` keeps it.
* `Panels ▾` (every plot tab) hides and shows the tab's panels; the panels
  left share the room. `Layers…` (Structure and Current tabs) hides parts of
  the 3-D view and sets the opacity of each part
  (see [BJT effects, breakdown and comparing designs](#bjt-effects-breakdown-and-comparing-designs)).
* `Current` → `Contact` adds the panel `Current under a contact`; the I-V tab
  has `Pin sweep`, `Pins…`, `A/B compare…` and `Table…` (the parameters of
  the last and the pinned sweeps side by side), and a BJT panel (gain and f_T
  against I_C) for Gummel sweeps.
* Region labels sit where they fit and never cover each other, the line cuts
  mark the regions they cross in a strip above the curves, and the semi-log
  I-V panel spans the currents that matter (a curve far below them is listed
  as "below the axis" instead of squeezing the plot).

**Settings and shortcuts**

* Text size: View → Text size, or Ctrl + / Ctrl − / Ctrl 0 (⌘ on a Mac). Plot
  labels scale too.
* Keys: F1 help, F5 run, F6 sweep, Esc stop, F11 fill the screen, Ctrl+E
  editor, Ctrl+O / Ctrl+S open / save the device, Ctrl+Q quit. On a Mac the
  Ctrl shortcuts use ⌘, and since F5 needs the fn key there, ⌘R runs and ⌘.
  stops.
* Help (F1) is a scrollable manual with a topic list and search.
* File → Save/Open device stores the whole editor state as JSON (defects and
  the adaptive mesh lines included). Export solution/sweep CSV and the
  visible plot (PNG/PDF/SVG) are also in the File menu.

**Gate stacks.** A gate on a semiconductor face takes an optional layer
stack in the Contacts page, semiconductor side first, in nm:
`SiO2 0.5, HfO2 2.5`. The gate capacitance then uses the stack's EOT, and
gate tunnelling sees every layer. A box gate (a metal electrode inside
meshed insulator regions) needs no stack.

**Nets.** Contacts with the same `Net` name are wired together outside the
simulated cell: the NMOS and PMOS gates of an inverter, the three sides of a
FinFET gate, a source tied to the substrate. A net carries one voltage, is
swept as one terminal, and its currents are summed in the plots.

**Floating terminal (I = 0).** In the Sweep page any contact or net can be
floating: at every point its voltage is solved so that its total current is
zero (a bracketed secant on the continued solution). This turns a sweep into
a transfer curve (the inverter's unloaded output). With `Adaptive
refinement`, steep parts of the curve get extra points automatically.

**Extracted parameters** (I-V and Results tabs): diode ideality and J₀ (a
straight line through the steepest-ideal part of ln J spanning two decades,
with kT/q at the sweep's temperature); for gate sweeps SS (the steepest
decade above the resolution floor), V_th (at |V_DS| ≤ 0.1 V the tangent at
maximum g_m minus V_DS/2, in saturation the √I_D extrapolation - Ortiz-Conde
et al. 2002 - plus constant current 100 nA·W/L where the template gives
W/L) and I_on/I_off; for gate-current
sweeps the peak J_G and J_G at ±0.5 and ±1 V; for transfer curves V_M, the
maximum gain, V_OH/V_OL, V_IL/V_IH, the noise margins and the peak supply
current; for output curves the Early voltage; for Gummel plots of a BJT the
current gain and f_T,max (f_T = ΔI_C / 2πΔQ between neighbouring points, at
their logarithmic-mean I_C); the avalanche breakdown voltage when impact
ionisation is on (a heuristic current knee otherwise). Pinned sweeps get the
same list.

**Open terminal.** A floating terminal whose own swept terminal is the one
observed (`Plot current of` = the swept terminal) is drawn as an I-V curve,
not a transfer curve, with the open terminal's voltage on the semi-log
panel's right axis: the open base of a BV_CEO measurement.

## Contact metals

Contacts page → `Metal`: Al, Ti, Ag, Cr, W, Mo, Cu, Au, Pd, Ni, Pt, TiN, TaN,
TiSi2, CoSi2, NiSi, PtSi, ErSi2, n⁺/p⁺ poly-Si, or `(ideal / custom φm)`.
The line under the field shows the barrier the metal forms on the
semiconductor under the contact and its source:

* **measured** barrier heights on n-type Si, GaAs, 4H-SiC and GaN (Sze & Ng;
  the SiC/GaN contact literature), e.g. Al/n-Si 0.72 eV, Ti/n-Si 0.50,
  PtSi/n-Si 0.84, Ni/4H-SiC 1.5 eV;
* the **Cowley-Sze pinning model** for any other pair:
  φBn = S(φm − χ) + (1 − S)φpin (Si S = 0.27, φpin 0.80 eV; GaAs strongly
  pinned, S ≈ 0.07; SiC and GaN S ≈ 0.4).

The hole barrier is Eg − φBn. On a semiconductor face both `Ohmic` and
`Schottky` contacts with a metal become physical metal contacts (Solver page,
two switches, on by default):

* **Thermionic emission** with the Richardson velocity v = A**T²/(qN_c,v)
  against the density in equilibrium with the metal, over a barrier lowered
  by the **image force** - flat up to x_m = Δφ/2F, −q/16πεx beyond - which is
  part of the band the carriers drift and diffuse in, screened by the free
  carriers (Debye length). Together this is Crowell-Sze
  thermionic-emission-diffusion; reverse currents keep rising with voltage as
  on real diodes.
* **Tunnelling through the barrier** (field and thermionic-field emission):
  Tsu-Esaki current with WKB transmission through the band along the contact
  normal, injected where the carriers leave the barrier. It makes a metal on
  heavily doped silicon ohmic: ρc ≈ 1e-7…1e-6 Ω·cm² on 1e20 cm⁻³.

A `Gate` with a metal takes the metal's work function. A Box contact keeps
an ideal ohmic contact (Ohmic) or the fixed barrier (Schottky). Real barriers
vary by ±0.1 eV with surface preparation: type φm to override.

## Defects

The `Defects` page puts defects into any device. A defect is a physical
effect from the library (literature values, every number editable) placed in
a shape.

**What a defect does**

* **Deep levels (traps):** SRH recombination through the level (any energy,
  its own capture cross sections, thermal velocities from the effective
  masses) and the charge it holds: an acceptor is negative when it holds an
  electron, a donor positive when empty. The occupancy follows the local
  electron and hole densities (steady state) and enters Poisson's equation -
  traps compensate doping, deplete regions and shift thresholds. A neutral
  level only recombines.
* **Interface traps, D_it** (cm⁻²eV⁻¹): a constant density across the gap at
  the semiconductor/insulator interfaces inside the shape and under face
  gates, as 24 levels (acceptor-like above midgap, donor-like below).
* **Fixed charge:** at interfaces (q/cm², oxide charge) or in a volume.
* **Inclusions:** the shape becomes another material - a void (ε = 1), SiO2
  (oxide precipitate) or a metal particle (ε → ∞, a floating equipotential) -
  with interface states on its surface. A metal particle's states pin the
  Fermi level at its Schottky barrier on that semiconductor: it depletes its
  surroundings like a floating Schottky contact and recombines like one.
* **Extended defects:** dislocations (states per cm of line), stacking
  faults and grain boundaries (states per cm², or a D_it band).
* **Radiation damage:** neutron fluence (the three-level Perugia model for
  p-type FZ Si, level densities η·Φ) and total ionizing dose (oxide charge +
  D_it).

**Library (26 presets)**

| Category | Presets |
|---|---|
| Voids & inclusions | void, oxide precipitate (SiO2), metal/silicide particle (any metal), Cu3Si precipitate, COP (crystal-originated void), micropipe (4H-SiC) |
| Point defects & contamination | custom deep trap, interstitial iron Fe_i (Ev + 0.38 eV), FeB pair (Ec − 0.23), gold (Ec − 0.55 / Ev + 0.35), platinum (Ec − 0.23 / Ev + 0.32) |
| Radiation damage | neutron fluence (Perugia), ionizing dose |
| Interfaces | Si/SiO2 D_it, 4H-SiC/SiO2 (D_it + near-interface acceptor traps), oxide fixed charge, volume charge |
| Extended defects | dislocation (Si), threading dislocation (GaN), stacking fault (4H-SiC), grain boundary (poly-Si) |
| SiC · GaN · GaAs | Z1/2, EH6/7 (SiC lifetime killers), C_N and Fe_Ga (semi-insulating GaN buffers), EL2 (GaAs) |

The description under the selection gives each value and its source.
Measured capture cross sections scatter by up to 10× (more for Pt): treat
them as starting points.

**Placement.** Shapes: sphere (diameter), box (dx, dy, dz), cylinder
(diameter, length along its axis, 0 = through the device), slab (thickness
along its normal axis, widths, 0 = across), everywhere (the whole host
material). `Host` limits a defect to one semiconductor. Manual: fill the form
and `+ Add`. Random: `Random scatter` places copies of the form's defect at
random centres in a box - a number, or a density (per cm³ of the box; per cm²
for lines through the device, e.g. a threading-dislocation density) - with a
seed (same seed, same layout) and a size spread; the copies form a group that
`Remove group` deletes at once.

**Exactness on any mesh.** The number of states is exact whatever the mesh:
deep levels in a shape smaller than the cells are spread over the nearest
nodes, a dislocation over the nearest row, and interface states carry the
true surface area of the shape (the cell faces around a meshed sphere add up
to ~1.5× its area). An inclusion replaces whole control volumes, so the
graded mesh is made denser through the first 40 inclusions, lines and planes
(boxes get exact faces); one smaller than the cells around it acts through
its surface states only, and the Defects page says so after a run.

**Background lifetime.** Solver page → `Background SRH τ0`: the lifetime of
the built-in midgap recombination (empty = material default, Si 10 µs; `1e-3`
= float-zone silicon; `none` = only the defects recombine).

**Seeing them.** The Structure tab draws every defect; solution maps outline
the ones the plane cuts (with `Regions` on); 3D Slicer quantities `trap
density` and `trapped charge`.

## Adaptive mesh refinement

Device page → `Refine the mesh to the solution on ▶ Run` (rectangular mesh).
After the solve, every node is checked along each axis: how far the
electrostatic potential, the band edges the carriers move in (with
band-gap narrowing and the quantum-confinement potential) and the
quasi-Fermi potentials of the carriers that matter there bend away from the
straight line through the two neighbours, in kT/q. Where that exceeds the
threshold (default 1 kT/q) the mesh does not resolve the solution - an
unresolved depletion edge, an inversion layer, current crowding - and the
interval is split into √(deviation/threshold) parts. The device is solved
again from the old solution interpolated onto the new mesh (from equilibrium
if that fails to converge), for up to `Passes` rounds or until the node
budget (`Max nodes ×` the unrefined count) is used; the largest deviations
are refined first. Uniform fields and resistors need no refinement and get
none. The added lines stay for sweeps and later runs until `Clear
refinement` or another template, and are saved with the device.

How it does (validate.py section 16, and a survey of the templates at their
default bias against meshes 2.5× finer in x and y):

| Device | Default mesh | After refinement (3 passes, ≤ 2× nodes) |
|---|---|---|
| p⁺/n diode on a coarse uniform mesh (15 nodes / 10 µm), I(0.5 V) | +3.0 % | +0.04 % (36 lines) |
| the same at 0.7 V | +3.6 % | +0.11 % |
| n-JFET template, drain current | +7.6 % | +0.9 % |
| NMOS template, linear region (V_G 1.5 V, V_D 50 mV) | +0.37 % | −0.06 % |
| PN Diode (+0.15 %) and Schottky (+0.5 %) templates | resolved | no lines added |
| PIN diode at 0.8 V (high injection) | −2.1 % | not flagged: no lines added |

It is not a cure-all: the subthreshold current of a short MOSFET with
quantum confinement depends on where the MLDA potential ends at the gate edge
(a model discontinuity no mesh resolves), and a few templates moved by a
similar small amount in the other direction (BJT −0.8 → +1.6 %). Compare
with a uniformly finer mesh when it matters.

## BJT effects, breakdown and comparing designs

**Impact ionisation (avalanche breakdown).** Solver page → `Impact
ionisation (avalanche breakdown)`. Carriers accelerated by a high field
create electron-hole pairs at G = (α_n|J_n| + α_p|J_p|)/q, with the
Chynoweth law α = γ a exp(−(γ b/F)^m), F the field along each carrier's
current and γ = tanh(ħω/2kT₀)/tanh(ħω/2kT) the optical-phonon temperature
factor (ħω = 63 meV, T₀ = 300 K), so breakdown voltages rise with
temperature. Coefficients: Si van Overstraeten - de Man (holes in two field
ranges split at 4e5 V/cm), GaAs Bulman et al., 4H-SiC Hatakeyama et al.
(along the c axis), GaN Ji et al.; Ge, InP and the alloys have none. The
generation is a source in both continuity equations, updated every Gummel
iteration, and enters the balance of floating regions. The 3D Slicer
quantity `log10 G avalanche` shows where the pairs are made.

Close to breakdown the multiplication grows without bound (~10-50) and the
Gummel iteration, which solves the equations one after the other, stops
converging. The I-V tab reports BV where 1/I, extrapolated from the last
converged points, reaches zero (1/I falls to zero at breakdown, for a
junction and for an open-base transistor alike) - or, with `BV knee` > 0,
where the current density reaches that value. Past breakdown the current is
held back only by the series resistance and the space charge of the
avalanche carriers; the solver may converge again on that high-current
branch (a jump in the curve). With impact ionisation on, a sweep ends after
two failed points in a row. Checked in validate.py section 17.

**Gain and f_T (BJT panel).** Every sweep point stores the electron and hole
charge in the device. For a Gummel sweep - the base swept, the collector
observed, V_CE held - the I-V tab's third panel draws β = I_C/I_B and the
quasi-static f_T = g_m / (2π C_in) against I_C, with g_m = dI_C/dV_BE and
C_in = dQ/dV_BE (Q: the charge of the base's majority carrier, holes for an
NPN): the charge-control estimate of where |h21| falls to 1. It contains the
emitter and collector depletion charges and the stored minority charge -
including the charge in the collector once the base pushes into it at high
current (the Kirk effect), which is what makes f_T fall past its peak. It
leaves out the non-quasi-static phase delay and every parasitic outside the
simulated cell (pads, interconnect, base resistance beyond the cell).

**Early voltage.** For output curves (collector or drain swept and observed,
three or more terminals) the flat part of the curve - after the saturation
knee, before any breakdown upturn - is fitted with a straight line: V_A =
I/g_o − V, where the dotted line meets I = 0 at −V_A (for a MOSFET, 1/λ).

**Current crowding.** Current tab → `Contact`: the panel `Current under a
contact` shows the current density through the chosen contact or net, node
by node - a profile across it in 2-D sections, a map over it in 3-D devices -
with the ratios J max/mean and centre/edge. Face contacts only.

**Comparing designs.**

* `Pin sweep` (I-V tab) keeps the last sweep as a reference curve. Later
  sweeps of the same swept terminal are drawn over it, dashed in the pin's
  colour, on the linear, semi-log and BJT panels, and its parameters are
  listed under those of the last sweep.
* `Pins…` shows or hides each pin (click `Shown`), renames it (double-click),
  changes its colour, removes it, exports all pins to CSV, and puts the
  fields of a pin's last bias point in the solution tabs (`Show fields`).
* `Table…` (I-V tab) lists the parameters of the shown pins and the last
  sweep side by side, one column each: the device, the sweep, the
  temperature and physics models (`default` = the Solver page's default
  switches, then `+` and `−` what differs), then every parameter the I-V tab
  extracts
  - V_th (max-g_m and constant-current), SS, g_m,max, I_on, I_off and their
  ratio; diode n and J_0; V_A and g_o; BV; β,max, f_T,max and the currents
  where they occur; V_M, V_OH/V_OL, V_IL/V_IH and the noise margins of an
  inverter - and the output current at the end of each sweep. Pick the
  `Reference` column; the other columns show their difference from it:
  voltages in mV or V (V_th, SS, V_M, n), currents and frequencies as ratios
  (×), gain, f_T, V_A and BV also in per cent. `Show` switches between values,
  values with differences and differences only; `Hidden pins too` adds the
  hidden pins. `Copy` puts the table on the clipboard tab-separated (paste it
  into a spreadsheet), `Export CSV…` writes it with the numbers in base
  units (V, A, S, Hz). A `*` marks a value with a note in the I-V tab (V_th
  or f_T still rising at the sweep end). The table follows every new sweep
  and pin.
* `A/B compare…` captures two designs and runs the Sweep page's sweep on both,
  one after the other. A design is the whole editor state - structure, mesh,
  contact voltages, physics switches, temperature and an open (floating)
  terminal - as it is when you press `Capture the editor's device`; `Load
  into the editor` brings it back. The templates of the `BJT effects &
  breakdown` group have variant pairs: choose one under `Variants`, press
  `→ A`, choose the other, `→ B`, then `▶ Run A and B`. Both results become
  pins (A blue, B pink); older pins are hidden. B can differ from A in
  anything: a doping, a finger layout, one physics switch, the temperature.

**Visibility.** `Panels ▾` on every plot tab hides panels so the others get
the room: the Structure tab's 3-D view, cross-section and legend; the 3D
Slicer's 3-D view and its three planes; the band diagram and ribbon; the five
Current panels; E-field and mobility; the three I-V panels. `Layers…` (on
the Structure and Current tabs) lists every region, contact and the defects
with its opacity. Click `Shown` to hide or show a part. Select one or more
parts (Ctrl or Shift for several) and move `Opacity of the selected parts` to
give them their own opacity, e.g. a see-through substrate under an opaque
channel; `Auto` hands them back to the sliders for all regions and all
contacts. `See-through contacts` sets the contacts to 0.15, `Insulators in
3-D` draws the oxides too, and `Apply to the cross-section too` uses the same
choices in the Structure tab's section. The Current tab's structure overlay
(`Map`, `Structure in 3-D`) follows the same list: hidden parts stay hidden,
the others take the tab's own `Opacity` unless they have one of their own.

**The demonstration templates** (group `BJT effects & breakdown (A/B
variants)`, 2-D sections, seconds per bias point):

| Template | Variants (A/B) | Sweep | What to look at |
|---|---|---|---|
| BJT effects: Kirk effect (f_T vs I_C) | collector 1e16 / 1e17 cm⁻³ | Base 0.6 → 1.1 V at V_CE = 2 V | f_T peaks at 6.3 GHz near 1.2e-5 A per cell (6.7e3 A/cm² over the emitter) with the light collector, 11.5 GHz near 1.6e-5 A with the heavy one; β 36 / 60 at low current |
| BJT effects: emitter current crowding | one 20 µm emitter / four 5 µm fingers | Base 0.6 → 0.98 V at V_CE = 2 V | Current tab, `Contact` = Emitter: at 0.98 V the middle of the wide emitter carries 1.4 % of its edge current density; the fingers deliver 3.5× its collector current for the same V_BE |
| BJT effects: Early effect (output curve) | base 5e16 / 2e17 cm⁻³ | Collector 0 → 5 V at V_BE = 0.75 V | V_A = 5.6 V (light base) / 42 V (heavy base, lower β) |
| Avalanche: Si p+n diode breakdown | impact ionisation on / off | Anode 0 → −56 V | BV ≈ 52 V (ionisation integral 51.9 V, Baliga 53.4 V); with impact ionisation off no breakdown at all |
| Avalanche: NPN BV_CBO, BV_CES, BV_CEO | base shorted to the emitter / base open / emitter left out (pick two) | Collector 0 → 60 V | BV_CBO ≈ 56 V > BV_CES ≈ 34 V (the avalanche holes flowing to the base contact turn the emitter edge on) > BV_CEO ≈ 21 V (open base: M = 1 + 1/β; Miller's BV_CBO/β^(1/4) gives 22 V) |

## Validation

`python3 validate.py` runs 158 checks against closed-form theory (about 10 min
on two cores); `--templates` adds section 9: all 213 checks in about 21 min.
Each check prints PASS/FAIL with its tolerance, and the script doubles as an
example of driving the core without the interface. The references are in
[SOURCES](SOURCES.md).

| § | Test | Result |
|---|---|---|
| 1 | Si p-n junction in equilibrium: built-in potential, np = nᵢ² | exact (1e-11 V, 1e-14) |
| 2 | Abrupt junction, 1–60 V: depletion width vs depletion approximation | within 0.05 % |
| 3 | Long-base p⁺/n diode vs Shockley theory, 0.45–0.55 V | within 0.5 % |
| 4 | 4H-SiC blocking junction to 2 kV: peak field | within 0.34 % |
| 5 | AlGaAs/GaAs and AlGaN/GaN heterojunctions: Vbi, ΔE_c = Δχ | exact |
| 6 | MOS capacitor: flat band, ψₛ(V_T) = 2φ_F, inversion charge | 0.06 mV; dQ/dV = 0.96 C_ox |
| 7 | Long-channel NMOS, linear region: I_D vs μQ_inv V_DS/L | +1.25 % |
| 8 | 2-D BJT (band-gap narrowing, τ(N)): base current vs Gummel number | within 0.45 % |
| 9 | (`--templates`) all 55 templates converge at their default bias | all converge |
| 10 | CMOS inverter, FinFET | V_T ±0.48 vs ±0.46 V; V_M 1.14 V; SS 67 mV/dec, DIBL 28 mV/V |
| 11 | Gate tunnelling: Tsu-Esaki kernel, thin SiO2, high-k, Fowler-Nordheim | kernel 0.01 %; 81 A/cm² at 1 V (1.2 nm); FN slope 2.25e10 V/m |
| 12 | MLDA quantum confinement, Lombardi mobility | ΔV_th +57 mV, centroid 0.82 nm; universal curve within 18 % |
| 13 | Prism mesh vs rectangular mesh; round MOS capacitor vs radial Poisson-Boltzmann | 1e-5 (all nodes), within 1 % (thinned); 1–4 mV |
| 14 | Schottky diodes vs Crowell-Sze; contact resistivity vs Tsu-Esaki; Al ohmic on n⁺, rectifying on n⁻ | within 0.22 %; within 0.3–5 % over 8 decades |
| 15 | Defects: compensation, trap-limited lifetime, D_it swing, oxide charge, void, metal particle, numbers of states | exact … 1.25 % |
| 16 | Adaptive refinement vs finer meshes; sweep predictor | coarse diode 3.0 % → 0.04 %; same currents, 671 vs 829 iterations |
| 17 | Avalanche breakdown of Si p⁺/n junctions vs the ionisation integral and Baliga; 400 K; stored charge; short-base diode transit time; Early voltage and f_T extraction | BV within 0.2 % of the ionisation integral (1e16, 1e17 cm⁻³, 400 K), 3 % under Baliga; BV ×1.21 at 400 K; charge sums exact; transit time −3 %; V_A, f_T, β within 0.03 % |
| 18 | Material parameters and analysis: Si mobility vs Masetti (1e15–1e20, 300 and 400 K), 4H-SiC holes vs Schaffer, nᵢ and V_bi at 250/400 K, the 32-contact limit, NaN from bad getter indices, ideality factor (n = 1 and 1.7, series resistance, 250–400 K), V_th (linear and saturation), SS above a noisy floor | mobility exact; n within 0.15 %; V_th exact; SS 0.0003 % |

## Device templates

The default bias is what **Run** solves; the sweep is what **Sweep** runs.
Currents are per simulated cell. Positive current flows into the device
through that terminal. y is depth: the top surface is y = 0 in every template.

| Template | What it shows | Default bias | Sweep |
|---|---|---|---|
| PN / PIN Diode, GaAs LED | Textbook diodes | 0 V | Anode −1…+1 V, forward positive |
| Schottky | W/n-Si diode (W contact metal: φB 0.67 eV) | 0 V | Metal −1.5…+0.7 V |
| **Metal contacts: ohmic vs Schottky (Al)** | Two aluminium contacts on n-Si: on an n⁺ implant (1e20) electrons tunnel through the 3 nm barrier - ohmic; on the n⁻ layer (2e16) only thermionic emission is left - a Schottky diode | 0 V | Right contact −2…+0.5 V: 4 nA at +0.3 V (ideality 1.04), −0.06 pA at −0.3 V, rectification 6e4 |
| NPN / PNP BJT | Vertical BJT with p⁺ extrinsic base | V_BE 0.7, V_CE 2 V | Base (Gummel plot of all terminals) |
| **Lateral PNP (bipolar IC)** | Junction-isolated lateral PNP (µA741-style): emitter and collector ring side by side in the n-epi, 2 µm lateral base, n⁺ buried layer against the parasitic vertical PNP | V_EB 0.7, collector −2 V, substrate −5 V | Emitter 0.3…0.8 V: β 40-50 falling to ~30 at high injection, substrate current ~2 % of I_C (~14 % without the buried layer) |
| **SOI lateral NPN** | Symmetric lateral NPN in a 60 nm SOI film (IBM-style CMOS-compatible bipolar), 100 nm base contacted through a p⁺ poly extrinsic base | V_CE 1 V | V_BE 0.5…1.0 V: ideal collector current, β 16-17 |
| NMOS / PMOS | 7 nm oxide MOSFETs, L = 320 nm, V_th ≈ 0.6 / −0.52 V; MLDA + Lombardi on | V_GS = V_DS = ±1.5 V | Drain 0…±3 V |
| **CMOS Inverter (0.25um)** | Twin-well NMOS + PMOS with STI; gates on net `Vin`, drains on net `Out` | Vin 0 V, Out floating | Vin 0…2.5 V with Out floating: transfer curve, gain ≈ 21, V_M ≈ 1.06 V |
| **FinFET (bulk tri-gate)** | 22 nm class: 10 nm fin, L_g 30 nm, 0.5 nm SiO2 + 2.5 nm HfO2, metal gate on three sides | V_G = V_D = 0.8 V | Gate 0…0.8 V: SS 67 mV/dec, I_on 48 µA/fin, I_on/I_off 3e6 |
| **FinFET (tapered, rounded fin)** | The same FinFET with a tapered fin and rounded corners, on the triangular mesh | V_G = V_D = 0.8 V | Gate 0…0.8 V |
| **GAA nanowire FET (round)** | 10 nm Si wire, gate all round, triangular mesh | V_G = V_D = 0.7 V | Gate 0…0.7 V: SS ≈ 65 mV/dec |
| Gate leakage: 1.2nm SiO2 / HfO2 high-k NMOS, Fowler-Nordheim: 8nm SiO2 MOS | Gate tunnelling demonstrators | V_G 1 V | Gate current sweeps |
| **Defects: leaky p-n junction (particle, void, dislocation)** | 3-D p⁺/n diode of float-zone quality with a NiSi precipitate in the depletion region, a void and a decorated dislocation | 0 V | Anode −2…+0.6 V: leakage −1.4e-14 A from the precipitate alone; I(0.3 V) ×1900 (particle), ×12 (void), ×7 (dislocation) |
| **Defects: NMOS after irradiation (oxide charge + D_it)** | The NMOS after a total ionizing dose: +5e11 q/cm² oxide charge, D_it 5e11 cm⁻²eV⁻¹ | V_G 1.5, V_D 0.05 V | Gate 0…1.5 V: V_th 0.574 → 0.501 V, SS 83.9 → 92.9 mV/dec, I_off ×18 |
| VDMOS, LDMOS | Power MOSFETs | V_GS 10 / 5 V, V_DS 5 / 15 V | Drain output curve |
| IGBT | Planar field-stop IGBT (600 V class) | V_GE 15, V_CE 2 V | Collector 0…5 V |
| GaAs HBT | Vertical AlGaAs/GaAs HBT | V_BE 1.35, V_CE 2 V | Collector output curve |
| Si Solar, 4H-SiC PiN | Dark I-V, kV blocking | 0 V | Forward / reverse |
| SiC JBS Diode | 4H-SiC JBS: Ti Schottky segments (1.1 eV) between ohmic p⁺ stripes | 0 V | −1200…+2 V |
| GaN HEMT | AlGaN/GaN polarisation 2DEG, Schottky gate | V_D 5 V | Drain 0…10 V |
| n-JFET | 0.8 µm channel between tied p⁺ gates | V_GS −0.5, V_DS 3 V | Drain 0…5 V |
| SiC / Si Trench MOS | Trench MOSFETs | on-state | Drain output curve |
| RW: 1N4148, 1N5819 (Mo barrier, φB 0.68 eV), 1N4733A, C4D | Rectifiers, Schottky, Zener, SiC JBS | 0 V | Reverse to rating … forward |
| RW: BPW34, Ge 1550 nm photodiodes | PIN photodiodes, dark | reverse bias | Anode sweep |
| RW: AlGaAs IR LED, PERC cell | DH emitter, thin-wafer PERC cell (dark) | 0 V | Forward positive |
| RW: BC547 | Small-signal NPN | V_BE 0.7, V_CE 3 V | Base |
| **RW: PNP / NPN power BJT 250 V (MJ21193 / MJ21194)** | Audio power transistors with a comb emitter: a 3-D section (0.28 × 0.20 mm) of interleaved emitter and base-contact fingers on bus bars, 7 µm base, 30 µm collector drift (5e14); the section plane opens as a top view of the two combs | V_BE ∓0.75, V_CE ∓5 V (PNP negative) | Base ∓0.45…−0.85 / +0.80 V: hFE 44 (PNP) / 56 (NPN), falling at high current as the collector goes into quasi-saturation (PNP 17 at ~17 A, NPN 46 at ~33 A for a 5 × 5 mm die, ×450 the section). Heavy: ~30 s per point |
| **RW: Lateral PNP LDO 1 A (LM2940)** | One round-emitter cell (simulated as a quarter, triangular mesh) of the lateral PNP pass transistor of a bipolar LDO: 12.7 µm p emitter, 5 µm lateral base in 3e15 n-epi, p collector ring, n⁺ buried layer, base contact through an n⁺ sinker | V_EB 0.85, V_EC 2 V, substrate −5 V: 2.8 mA per emitter at β 25 | Collector: V_EC 2 → 0.05 V (dropout). Below ~0.2 V the transistor saturates: I_C −60 % at 0.05 V, base current ×4, substrate current ×10 |
| RW: IRLZ44N, superjunction 600 V, SiC 900 V (C3M) | Power MOSFETs in the on-state | V_DS 1–2 V | Drain output curve |
| RW: IGBT 1200 V (IKW40N120) | Field-stop trench IGBT with a transparent collector | V_GE 15, V_CE 1.5 V | Collector 0…1.8 V |
| RW: GaN RF HEMT 28 V, GaN HEMT 650 V (GS66508) | RF HEMT; e-mode p-GaN gate | V_D 28 / 1 V | Drain / gate |
| RW: Thyristor 800 V (BT151) | 4-layer pnpn; latching near V_G ≈ 0.65 V | V_AK 2 V | Gate 0…1 V |

For the blocking state of a power switch, set its gate to 0 V and sweep the
drain or collector to the rated voltage. The field distribution and E_max are
meaningful; switch on impact ionisation (Solver page) for the breakdown
voltage itself.

The group `BJT effects & breakdown (A/B variants)` - Kirk effect, emitter
current crowding, Early effect, Si avalanche diode, NPN BV_CBO / BV_CES / BV_CEO - is
described in [BJT effects, breakdown and comparing designs](#bjt-effects-breakdown-and-comparing-designs).

## How the core works

* **Meshes.** Every equation is assembled with the box method: 7-point
  stencils on the rectangular mesh, a general CSR graph (edge lengths, dual
  face areas, node volumes) on a prism mesh. The same Scharfetter-Gummel,
  Poisson, MLDA, Lombardi, tunnelling, metal-contact and trap code runs on
  both.
* **Poisson.** Damped Newton with line search; CG with a symmetric
  aggregation-multigrid preconditioner. Trapped charge and its derivative
  enter the Newton step.
* **Continuity.** Scharfetter–Gummel fluxes with band-modified potentials
  (heterojunctions, band-gap narrowing, quantum confinement, image-force
  band), solved in scaled form so densities spanning 40+ decades keep uniform
  relative accuracy; BiCGSTAB with a Petrov–Galerkin aggregation multigrid.
* **Coupling.** A Gummel map on (ψ, φₙ, φₚ) with Anderson acceleration
  (depth 12).
* **Continuation.** Adaptive voltage steps; inside a ramp each step starts
  from the extrapolation of the last two. Steps that converge slowly (more
  than 25 iterations) are followed by 2.2× longer ones: there the iteration
  count is set by the Gummel contraction rate, not by the step length, so
  fewer, longer steps are cheaper (−24 % iterations over the default solves
  of the 47 templates of v0.1; e.g. the GaN RF HEMT 2.2× and the Si trench MOSFET 3.7×
  faster). MOS gates are ramped before the drain or collector.
* **Sweep predictor.** When a solve continues the previous bias change in the
  same direction - the next point of a sweep, or the secant steps of a
  floating terminal - it starts from the secant extrapolation (or
  interpolation) of the states before and after the previous solve instead of
  from the old state: an O(ΔV²) instead of O(ΔV) start, 20–30 % fewer Gummel
  iterations on typical sweeps and about half on the IGBT output curve, with
  the same results to ~1e-5.
* **Metal contacts.** Per contact node, the thermionic-emission exchange with
  the metal (Richardson velocity against the density in equilibrium with
  the metal) is implicit in the continuity matrix. The image-force barrier is
  a band modification (flat to x_m, then −q/16πεx, screened by the local
  Debye length) on the carrier pushed away by the field. Tunnelling: per
  contact-node column (60 nm along the normal), the Tsu-Esaki current with
  WKB transmission through the actual band, deposited at the classical
  turning point (linearly between the two nodes around it). A velocity-
  saturation taper switches saturation off for motion across the barrier.
* **Deep levels.** Every trap species is a level (type, energy, σₙ, σₚ); a
  node carries any number of species with their densities. SRH through each
  level with n₁ and p₁ from its energy (relative to Ec, Ev or midgap),
  occupancy f = (cₙn + cₚp₁)/(cₙ(n+n₁) + cₚ(p+p₁)), charge −N_t f (acceptor)
  or +N_t (1−f) (donor), and the derivatives for Poisson's Newton step and the
  continuity matrices. Ohmic contacts use the neutral potential including
  the trapped charge.
* **Adaptive mesh refinement** (Python, `amr_run`): the indicator described
  above, new mesh lines, the old solution moved onto the new mesh by
  tri-linear interpolation (`api3_get_state` / `api3_set_state`), and a
  resolve without a bias ramp.
* **Impact ionisation.** At every semiconductor node G = (α_n|J_n| +
  α_p|J_p|)/q with the node currents from the face fluxes and F = |E·J|/|J|
  (the field along each current), updated every Gummel iteration as an
  explicit source of both continuity equations and of the floating-region
  balance; Chynoweth α with the optical-phonon temperature factor.
* **Stored charge.** `api3_carrier_sums` returns Σn dV and Σp dV; every sweep
  point records them (the quasi-static f_T).
* **Maximum principle, floating regions, mobility, recombination, MLDA,
  Lombardi, gate tunnelling, stop/progress, terminal currents:**
  quasi-Fermi potentials clamped to the span of the terminal voltages, a
  balance step for floating quasi-neutral regions, Masetti (Si) and
  Caughey-Thomas (other materials) mobility with velocity saturation, SRH τ(N) + Auger + radiative, Klaassen BGN, MLDA
  per-valley DOS cut at interfaces, Lombardi surface mobility, Tsu-Esaki/WKB
  gate tunnelling, cooperative stop, Kirchhoff-exact terminal currents with
  a round-off resolution floor.

## Limitations

* **Impact ionisation is local** (the field at each node - no dead space, no
  energy relaxation): high-field regions shorter than ~100 nm break down at
  too low a voltage. Near breakdown the Gummel iteration fails, so the
  breakdown voltage comes from the current's runaway; snapback, second
  breakdown and avalanche ruggedness are out of reach. The runaway shows only
  when the leakage it multiplies lies above the current resolution: in very
  small cells with long lifetimes, shorten the background lifetime (Solver
  page, τ0) - it raises the leakage, not the breakdown voltage. No
  coefficients for Ge, InP and the alloys. Off by default (Solver page).
* **No band-to-band tunnelling.** Zener breakdown and GIDL are absent. It
  also adds to avalanche breakdown above ~5e16 cm⁻³, and the silicon
  ionisation coefficients are fitted for 1.75–6e5 V/cm: an abrupt Si junction
  at 1e17 cm⁻³ breaks down at 14 V here against about 11 V in the classical
  curves (Sze & Gibbons); at 1e16 cm⁻³ the 52 V agree with the literature.
* **Currents near avalanche** are converged to their resolution floor; when
  the multiplied current is below ~1e-16 A per cell it can still move by a
  few per cent between solves (the floor in the results says by how much).
* **Boltzmann statistics, full ionisation.** No Fermi-Dirac statistics
  (degenerate n++/p++ layers, GaAs above ~4e17 cm⁻³ and HEMT channels are
  approximate) and every dopant is ionised (p-type 4H-SiC with Al and p-GaN
  with Mg are far too conductive).
* **No negative differential mobility.** GaAs and InP electrons saturate
  smoothly at v_sat (no transferred-electron peak, no Gunn effect): short
  GaAs/InP channels carry 2–3× too little current near the peak field.
* **Temperature.** Silicon's saturation velocity follows Canali et al.; the
  other materials keep their 300 K saturation velocity, and the SiC and GaN
  Auger and ionisation coefficients do not depend on temperature. Silicon's
  nᵢ is 1.07e10 cm⁻³ at 300 K (Green's densities of states), 11 % above the
  modern 9.65e9: saturation currents are ~23 % higher.
* **f_T is quasi-static** (charge control, no AC analysis): no
  non-quasi-static phase delay, no parasitics outside the simulated cell.
* **No optical generation.** Solar cells and photodiodes are simulated dark.
* **Quantum confinement is MLDA only** (no subband energies); heterostructure
  2DEGs stay classical. MLDA ends abruptly where a face gate ends, which makes
  subthreshold currents of short MOSFETs mesh-sensitive (refine uniformly to
  check; adaptive refinement does not fix a model discontinuity).
* **Gate tunnelling:** no image-force lowering in the oxide, no trap-assisted
  tunnelling, no valence-electron tunnelling.
* **Metal contacts:** one barrier height per metal/semiconductor pair (real
  barriers vary by ±0.1 eV with preparation - override with φm); no
  interfacial layer, no thermionic emission at heterointerfaces (an abrupt
  HBT therefore has a low β).
* **Defects:** traps are in steady state - no transient trapping, current
  collapse or hysteresis; no field-enhanced emission (Poole-Frenkel,
  trap-assisted tunnelling); a metal particle carries no current through
  itself (one spanning a junction does not short it); dislocations are not
  conducting pipes. Capture cross sections are literature values with large
  scatter.
* **Drift-diffusion transport.** No ballistic transport, strain or
  remote-phonon scattering; Lombardi has Si parameters only.
* **Not modelled:** self-heating, field plates, transients/AC.
* **Adaptive mesh refinement** works on the rectangular mesh only; a new line
  runs through the whole device. It resolves what the indicator sees
  (space charge, band bending, quasi-Fermi curvature); quantities that hinge
  on something else are not improved.
* **Low-current accuracy.** Currents are taken at the contacts, where they
  are a small difference of large drift and diffusion terms: below about
  1e-2 A/cm² they are good to 1e-4…1e-3 (they change by that much with the
  bias path or the number of threads), near the resolution floor to about
  0.1–1 %, and below it they are noise (the GUI marks them).
* **Limits.** At most 32 contacts and 16 different materials (gate-stack and
  defect materials count); more gives an error.
* **Speed.** Gummel iteration is slow under strong high-level injection: the
  1200 V IGBT takes about 2 min at its default bias on two cores, most
  templates seconds. The Gummel loop stops converging in a few high-current
  regimes (the 1200 V IGBT near V_CE ≈ 2 V; trench MOSFETs at full gate drive
  beyond V_DS ≈ 3–7 V; the superjunction beyond V_DS ≈ 5 V); the default
  sweeps stay inside the converging range.
* **Memory.** About 760 bytes per node (a prism mesh about 30 % more). The GUI
  refuses meshes larger than 80 % of the free RAM.
* **Triangular mesh.** A 2-D triangulation extruded along one axis, not a
  general tetrahedral mesh.
