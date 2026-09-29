# SemiSim 3D v7.4

A 3-D drift-diffusion simulator for semiconductor devices: a C core
(`semiconductor_3d.c` → `semiconductor_core.so`) driven by a Python/Tk GUI
(`main.py`). Devices are built from material/doping blocks - boxes, or boxes
cut to a round or sloped outline - on one of two meshes:

* **Rectangular:** a non-uniform tensor grid, refined automatically at
  junctions, gates, heterointerfaces and sheet charges.
* **Triangular (prisms):** a Delaunay triangulation of the cross-section,
  extruded along one axis (`semimesh.py`). It meshes round and sloped
  interfaces as they are (nanowires, rounded and tapered fins) and keeps fine
  rows only where the device needs them.

The core solves the coupled Poisson and electron/hole continuity equations
with the box method on either mesh and returns bands, quasi-Fermi levels,
fields, current densities, mobilities, recombination and terminal currents.

`validate.py` checks the core against closed-form semiconductor theory. The
42 built-in templates range from textbook structures, through a CMOS inverter,
bulk FinFETs (sharp, and tapered with rounded corners), a gate-all-around
nanowire FET and gate-leakage demonstrators, to representative commercial
parts. They are built from published device-class design rules; they are not
reverse-engineered dies.

For nanoscale MOS devices the core adds quantum confinement at oxide
interfaces (MLDA), surface mobility (Lombardi) and gate tunnelling through
any stack of SiO2, HfO2, Al2O3 and Si3N4 (Tsu-Esaki with WKB transmission:
direct and Fowler-Nordheim tunnelling).

## Files

| File | Purpose |
|---|---|
| `main.py` | GUI, device templates, rectangular mesh generator, experiment engine (nets, sweeps, floating nodes), ctypes bindings |
| `semimesh.py` | triangular-prism mesher: outlines (shapes), point placement, Delaunay box-method geometry, contacts, interface walls, display resampling. Must sit next to `main.py`. |
| `semiconductor_3d.c` | simulation core (single file) |
| `semiconductor_core.so` | portable Linux x86-64 build (x86-64-v2 baseline: any 64-bit Intel/AMD CPU, also under WSL2). `build.sh` rebuilds it for your CPU and runs about 15 % faster. |
| `validate.py` | physics validation battery (sections 1–13, PASS/FAIL) |
| `build.sh` | one-line build of the core |

## Build and run

Requirements: a C compiler with OpenMP, and Python 3.9+ with `numpy`,
`matplotlib` and `tkinter` (on Debian/Ubuntu the last one is `python3-tk`).

The shipped `semiconductor_core.so` runs as-is on x86-64 Linux, WSL2 included:
just `python3 main.py`. To build for your own CPU (Linux, gcc):

```bash
gcc -O3 -march=native -ffast-math -fno-finite-math-only -fopenmp -shared -fPIC \
    -o semiconductor_core.so semiconductor_3d.c -lm        # or: ./build.sh
python3 main.py
```

Do not copy a `-march=native` build to a different machine: it can die with
"Illegal instruction", because the sandbox CPU used here has AVX-512 and most
laptop CPUs do not.

* **macOS:** use Homebrew `libomp`:
  `clang -O3 -ffast-math -fno-finite-math-only -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -shared -fPIC -o semiconductor_core.so semiconductor_3d.c`
* **Windows:** build under MSYS2 MinGW-w64:
  `gcc ... -shared -o semiconductor_core.dll semiconductor_3d.c`.
  Then set `SEMISIM_SO=path\to\semiconductor_core.dll` before starting `main.py`.
* **`-fno-finite-math-only` is required.** Without it, `-ffast-math` deletes
  the solver's NaN guards.
* **Threads:** the GUI uses every core (Solver page → CPU threads).
  `OMP_NUM_THREADS` sets the default. A short OpenMP spin
  (`GOMP_SPINCOUNT=30000`) is preset, so other busy programs slow the solver
  down less.

## Using the GUI

The window adapts to the screen: 1280×720 and larger work, and the layout is
remembered between sessions (in `~/.semisim3d.json`).

```
┌ File  Devices  Run  View  Help ──────────────────────────────────────────┐
│ ☰  Device [CMOS Inverter ▾]  ▶ Run  ↔ Sweep  ■ Stop  [████░░] point 7/26… ? │
├──────────── editor (drawer) ──┬── plots ──────────────────────────────────┤
│ Device│Regions│Contacts│Solver│ Structure│3D Slicer│Bands│Current│I-V│…   │
│ Sweep                         │                                           │
│  (every page scrolls both     │   (matplotlib toolbar under every plot)   │
│   ways; wheel / Shift+wheel)  │                                           │
└───────────────────────────────┴───────────────────────────────────────────┘
```

**Toolbar**

* `☰` (Ctrl+E) hides the editor so the plots get the whole window. Drag the
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
  (`Rectangular` or `Triangular (prisms)`, with the prism axis) and the final
  mesh statistics (nodes, memory, minimum spacing or triangle quality).
* `Regions` and `Contacts`: scrollable tables. Select a row to edit it.
* `Solver`: temperature, physical models (band-gap narrowing, τ(N),
  velocity saturation, quantum confinement, Lombardi surface mobility, gate
  tunnelling), iteration limits, tolerances and CPU threads. Each template
  sets the model switches it was designed with.
* `Sweep`: the swept terminal, range and points, a floating terminal and the
  plotted current.

**Plot tabs**

* Structure (3-D view plus an XY/YZ/XZ cross-section with the real mesh -
  grid lines, or the triangles of a prism mesh; outlines drawn as they are;
  hover to read the region under the cursor), 3-D Slicer, Bands, Current,
  I-V, and Field & Mobility.
* `Results`: every terminal's applied and solved voltage, current, current
  density and resolution, the parameters extracted from the last sweep, and a
  run log.

**Settings and shortcuts**

* Text size: View → Text size, or Ctrl + / Ctrl − / Ctrl 0. Plot labels
  scale too.
* Help (F1) is a scrollable manual with a topic list and search.
* File → Save/Open device stores the whole editor state as JSON. Export
  solution/sweep CSV and the visible plot (PNG/PDF/SVG) are also in the File
  menu.

**Gate stacks.** A gate on a semiconductor face takes an optional layer
stack in the Contacts page, semiconductor side first, in nm:
`SiO2 0.5, HfO2 2.5`. The gate capacitance then uses the stack's EOT, and
gate tunnelling sees every layer. A box gate (a metal electrode inside
meshed insulator regions) needs no stack: tunnelling follows the meshed
layers between the semiconductor and the metal, as in the FinFET.

**Nets.** Contacts with the same `Net` name are wired together outside the
simulated cell, for example:

* the NMOS and PMOS gates of an inverter;
* the three sides of a FinFET gate;
* a source tied to the substrate.

A net always carries one voltage, is swept as one terminal, and its currents
are summed in the plots. Editing one member's voltage updates the whole net.

**Floating terminal (I = 0).** In the Sweep page any contact or net can be
floating. At every point its voltage is solved so that its total current is
zero, using a bracketed secant on the continued solution. The bracket comes
from the maximum principle: the root lies between the lowest and highest
applied voltage. This turns a sweep into a transfer curve (the inverter's
unloaded output). With `Adaptive refinement`, steep parts of the curve get
extra points automatically.

**Extracted parameters** (I-V and Results tabs):

* diode ideality factor and J₀;
* for gate sweeps: subthreshold swing SS, V_th (maximum-g_m extrapolation,
  plus constant current 100 nA·W/L where the template gives W/L) and
  I_on/I_off;
* for gate-current sweeps: the peak gate current and J_G at ±0.5 and ±1 V;
* for transfer curves: V_M, the maximum gain, V_OH/V_OL, V_IL/V_IH, the noise
  margins and the peak supply current;
* a heuristic current knee for blocking sweeps.

## Triangular (prism) mesh

Device page → `Mesh: Triangular (prisms)`. The cross-section normal to the
`Prism axis` is triangulated (Delaunay) and extruded along the axis over the
same 1-D grid the rectangular mesh uses there. The control volumes are the
Voronoi cells of the triangles times the 1-D cells, so every edge couples its
nodes through (dual face area / edge length) exactly as on the rectangular
mesh, and the physics is unchanged. `auto` picks the axis of the template's
outlines, else the axis with the fewest base nodes (z for a 2-D section).

**Point placement.**

* *Box geometry.* Near every interface, junction, gate and contact edge the
  points are the rectangular mesh's own nodes (its interface pairs and
  inversion-layer rows included), so both meshes resolve interfaces the same
  way. Away from them a point survives only if the local target spacing
  (growing with the distance to the nearest feature) still needs its lines:
  the fine lines a tensor grid drags through the whole device are dropped
  there. Typically 10–50 % fewer nodes; results agree with the rectangular
  mesh within ~1 %.
* *Outlines.* Rows of points along the outline's normals: a node pair
  straddling each interface (0.15 nm from a semiconductor interface, so the
  face between the two nodes IS the interface), a geometric boundary layer
  into the semiconductor (0.5 nm first step, ×1.3), the metal row exactly on
  a wrapped gate's surface. All members of an offset family (a gate stack)
  share the normals, so the layers stay aligned node by node and tunnelling
  paths run straight through them.

**Outlines (shapes).** A region or a Box contact can carry a
`semimesh.Shape`: a convex polygon with rounded corners, or a circle, in the
plane normal to one axis. The region is its box cut to the outline; a gate
electrode is its box minus the outline (`hole`). `shape.offset(t)` gives the
parallel outline at distance t (a conformal film: corner radii r + t). The
new templates use them; the tables mark them ◯, and Save/Open keeps them. On
the rectangular mesh an outline becomes a staircase of nodes.

**Display.** Line cuts and maps use the rectangular lines as a display grid:
every prism node is on it or interpolated inside its triangle (carrier
quantities from semiconductor nodes only, band edges from the node's own
material). In the 3D Slicer the map of the cross-section normal to the prism
axis is drawn on the triangulation itself. File → Export solution CSV writes
the prism nodes themselves.

**Speed.** The prism graph has its own multigrid (strength-based pairwise
aggregation). It is less effective than the rectangular one, so box geometry
runs at a similar speed with fewer nodes - faster on some large power
devices (the 1200 V IGBT about 5×), slower on others (VDMOS, SiC MOSFET about
1.5–2×). Use it where the geometry needs it.

## Validation (`python3 validate.py`)

All 78 checks of the default run pass (about 6 min on two cores);
`--templates` adds section 9 (every template converges, 42 more checks:
120/120 in about 13 min).

| § | Test | Result |
|---|---|---|
| 1 | Si p-n junction in equilibrium: built-in potential; mass action np = nᵢ² | Vbi exact to 1e-11 V; max\|np/nᵢ²−1\| = 1e-14 |
| 2 | Reverse-biased abrupt junction, 1–60 V: depletion width vs depletion approximation | −0.05 … −0.001 % |
| 3 | Long-base p⁺/n diode vs Shockley theory (coth form), 0.45–0.55 V | +0.45 / +0.08 / −0.19 %; at 0.35–0.40 V the +1…3 % excess is depletion-region SRH current, as expected |
| 4 | 4H-SiC p⁺/n⁻/n⁺ blocking junction to 2 kV: peak field | −0.34 / −0.11 / −0.08 % at 100 / 1000 / 2000 V |
| 5 | AlGaAs/GaAs and AlGaN/GaN heterojunctions: Vbi and ΔE_c = Δχ | exact |
| 6 | MOS capacitor (p-Si 1e17, 10 nm oxide): ψₛ at V_FB, ψₛ(V_T) = 2φ_F, dQ_inv/dV_G | 2e-16 V; −0.06 mV; 0.961·C_ox (inversion-layer capacitance) |
| 7 | Long-channel NMOS, linear region: I_D vs μQ_inv V_DS/L | +1.25 %; Kirchhoff 5e-10 |
| 8 | 2-D n⁺pn BJT with band-gap narrowing and τ(N): base electron current vs Gummel-number integral | +0.09 … +0.44 %; Kirchhoff ≤ 3e-7 (≥ 0.55 V); β = 43 |
| 9 | (`--templates`) every GUI template converges at its default bias, with its own model switches and mesh | all 42 converge |
| 10 | CMOS inverter and FinFET (12 checks) | see below |
| 11 | Gate tunnelling (9 checks) | see below |
| 12 | Quantum confinement (MLDA) and Lombardi mobility (5 checks) | see below |
| 13 | Triangular-prism mesh (19 checks) | see below |

Section 10 details:

* **CMOS inverter, thresholds** (classical models, against long-channel
  theory): V_th +0.484 / −0.487 V vs the estimate ±0.460 V.
* **CMOS inverter, transfer curve:** every point solved with the output
  floating; V_OH = 2.5 V, V_OL ≈ 1 nV. The output current is balanced to
  0.1 % of the supply current. V_M = 1.136 V (square-law estimate 1.07 V),
  and Vout is monotonic.
* **FinFET** (as the GUI runs it: MLDA, Lombardi, tunnelling, SiO2/HfO2
  stack): SS = 67.0 mV/dec (the thermionic limit is 59.6), DIBL 28 mV/V,
  I_on/I_off = 3.1e6, gate leakage 2.8 pA at V_G = V_D = 0.8 V.

Section 11 details:

* **Tunnelling kernel** against an independent integration (numerical WKB
  in x, fine energy grid): direct tunnelling through 1.2 nm SiO2, both
  directions, and a Fowler-Nordheim (triangular) barrier: all within 0.01 %.
* **n-MOSFET, 1.2 nm SiO2, V_G = 1 V:** J_G = 81 A/cm² (textbook range
  1e2–1e3 A/cm²); 5.0 decades per nm of oxide between 1.0 and 1.5 nm (the
  classic "10× per 0.2 nm"); the same J_G with and without MLDA (ratio
  0.99); Kirchhoff with the gate current closed to the round-off floor.
* **High-k:** 0.5 nm SiO2 + 3.95 nm HfO2 (same EOT) leaks 5700× less.
* **Fowler-Nordheim** through 8 nm SiO2, 7–12 MV/cm: the slope of the
  ln(J/E²) vs 1/E plot is B = 2.25e10 V/m, in the measured Si/SiO2 range
  (2.3–2.5e10; the barrier value for 3.1 eV and m_ox = 0.42 is 2.42e10 - the
  degenerate inversion-layer supply lowers the fitted slope slightly).

Section 12 details:

* **MLDA** wall factor against direct integration (both Si valleys): 2e-11.
* **MOS capacitor, 2 nm SiO2, p 5e17:** quantum threshold shift +57 mV
  (Schrödinger-Poisson literature: 40–100 mV); inversion-charge centroid at
  1e13 cm⁻² 0.82 nm (classical 0.37 nm); the capacitance-equivalent thickness
  grows by 0.08 nm (full Schrödinger-Poisson gives 0.1–0.3 nm: MLDA is a
  near-wall correction and captures the lower end).
* **Lombardi:** long-channel NMOS (5 nm SiO2, p 1e17), μ_eff against the
  universal curve 540/(1+(E_eff/0.9 MV/cm)^1.85) at 0.47–0.99 MV/cm: within
  6–18 %.

Section 13 details:

* **Prism mesh with every rectangular node kept** (the Delaunay triangulation
  of the tensor grid, solved through the general graph path): the 1.2 nm
  SiO2 NMOS (MLDA, tunnelling) and the FinFET (3-D, box gates) reproduce the
  rectangular solution to better than 1e-5.
* **Thinned prism mesh against the rectangular mesh** (same physics, each
  template's own models): NMOS I_D −0.31 %, 1.2 nm SiO2 NMOS gate current
  −0.003 %, GaN HEMT I_D −0.007 %, FinFET I_D −0.71 %, with 18–24 % fewer
  nodes.
* **Round MOS capacitor** (Si wire, R = 30 nm, p 1e19, in 2 nm SiO2 inside a
  wrapped metal gate; body contact on the axis) against the radial
  Poisson-Boltzmann solution: surface band bending within 1.1 / 1.2 / 2.1 /
  0.4 / 4.0 mV from accumulation to strong inversion (V_G −1.5 … +2.5 V). A
  rectangular (staircase) mesh of similar size misses by 3 / 19 / 75 / 54 /
  77 mV - the round mesh is 19× more accurate.
* **GAA nanowire FET** (10 nm wire, L_g 20 nm): SS 65 mV/dec, I_on/I_off
  3.6e6, I_on 11 µA at 0.7 V, gate leakage 2 pA.

The tolerances are in the script. It builds each device directly on the core
(sections 1–8, 11–12) or through the same engine the GUI uses (sections 10
and 13),
so it doubles as example code for scripting the simulator without the GUI.

## Device templates

The default bias is what **Run** solves; the sweep is what **Sweep** runs.
Currents are per simulated cell. Positive current flows into the device
through that terminal. y is depth: the top surface is y = 0 in every template.

| Template | What it shows | Default bias | Sweep |
|---|---|---|---|
| PN / PIN Diode, Schottky, GaAs LED | Textbook diodes | 0 V | Anode (or metal) −1…+1 V, forward positive |
| NPN / PNP BJT | Vertical BJT with p⁺ extrinsic base | V_BE 0.7, V_CE 2 V | Base (Gummel plot of all terminals) |
| NMOS / PMOS | 7 nm oxide MOSFETs, L = 320 nm, V_th ≈ ±0.6 / −0.52 V; MLDA + Lombardi on | V_GS = V_DS = ±1.5 V | Drain 0…±3 V |
| **CMOS Inverter (0.25um)** | Twin-well NMOS + PMOS with STI; gates on net `Vin`, drains on net `Out`; MLDA + Lombardi on | Vin 0 V, Out floating | Vin 0…2.5 V with Out floating: transfer curve, gain ≈ 21, V_M ≈ 1.06 V, noise margins, supply-current peak ≈ 15 µA |
| **FinFET (bulk tri-gate)** | 22 nm class: 10 nm fin, L_g 30 nm, 0.5 nm SiO2 + 2.5 nm HfO2 (EOT 0.94 nm), metal gate on three sides (net `Gate`), punch-through stopper; MLDA, Lombardi, tunnelling | V_G = V_D = 0.8 V | Gate 0…0.8 V: SS 67 mV/dec, V_th 0.29 V, I_on 48 µA/fin, I_on/I_off 3e6, I_G 3 pA |
| **FinFET (tapered, rounded fin)** | Same FinFET with a realistic fin: 8 nm at the top, 12 nm at the STI (86.7° sidewalls), 3 nm corner radius, conformal SiO2/HfO2 stack, wrapped metal gate. Opens on the triangular mesh. | V_G = V_D = 0.8 V | Gate 0…0.8 V: I_on ≈ 45 µA/fin |
| **GAA nanowire FET (round)** | 10 nm Si wire, L_g 20 nm, gate all round: 0.6 nm SiO2 + 1.9 nm HfO2, metal 4.55 eV; MLDA, Lombardi, tunnelling. Opens on the triangular mesh. | V_G = V_D = 0.7 V | Gate 0…0.7 V: SS ≈ 65 mV/dec, I_on ≈ 11 µA, I_on/I_off ≈ 4e6 |
| **Gate leakage: 1.2nm SiO2 NMOS** | Direct tunnelling through a 1.2 nm SiO2 gate (metal gate 4.15 eV, p 1e18) | V_G 1 V | Gate −1.5…+1.5 V, plotted current = gate: J_G(+1 V) ≈ 80 A/cm² |
| **Gate leakage: HfO2 high-k NMOS** | Same transistor and EOT with 0.5 nm SiO2 + 3.95 nm HfO2 | V_G 1 V | Gate −1.5…+1.5 V: J_G(+1 V) ≈ 1.5e-2 A/cm² (~5000× less) |
| **Fowler-Nordheim: 8nm SiO2 MOS** | Thick-oxide Fowler-Nordheim tunnelling | V_G 1 V | Gate 4…10 V: J_G 1e-9 → 0.4 A/cm² |
| VDMOS, LDMOS | Power MOSFETs | V_GS 10 / 5 V, V_DS 5 / 15 V | Drain output curve |
| IGBT | Planar field-stop IGBT (600 V class): JFET implant, conductivity modulation | V_GE 15, V_CE 2 V | Collector 0…5 V (0.7 V knee) |
| GaAs HBT | Vertical AlGaAs/GaAs HBT: 5e17 emitter over a 4e19 base | V_BE 1.35, V_CE 2 V | Collector output curve |
| Si Solar, 4H-SiC PiN, SiC JBS | Dark I-V, kV blocking | 0 V | Forward / reverse |
| GaN HEMT | AlGaN/GaN polarisation 2DEG, Schottky gate | V_D 5 V | Drain 0…10 V |
| n-JFET | 0.8 µm channel between tied p⁺ gates, V_GS(off) ≈ −1.6 V | V_GS −0.5, V_DS 3 V | Drain 0…5 V |
| SiC / Si Trench MOS | Trench MOSFETs (SiC: shield under the trench bottom, tied to source; current-spreading layer) | on-state, V_DS 2 / 1 V | Drain 0…5 / 0…3 V |
| RW: 1N4148, 1N5819, 1N4733A, C4D | Rectifiers, Zener, SiC JBS | 0 V | Reverse to rating … forward |
| RW: BPW34, Ge 1550 nm photodiodes | PIN photodiodes, dark (no optical generation) | reverse bias | Anode sweep |
| RW: AlGaAs IR LED, PERC cell | DH emitter, thin-wafer PERC cell (dark) | 0 V | Forward positive |
| RW: BC547 | Small-signal NPN | V_BE 0.7, V_CE 3 V | Base |
| RW: IRLZ44N, superjunction 600 V, SiC 900 V (C3M) | Power MOSFETs in the on-state | V_DS 1–2 V | Drain output curve |
| RW: IGBT 1200 V (IKW40N120) | Field-stop trench IGBT with a transparent collector: about 235 A/cm² at 1.5 V, 2–3× a real part (no channel-mobility degradation). Heavy: about 6 min on one core. | V_GE 15, V_CE 1.5 V | Collector 0…1.8 V |
| RW: GaN RF HEMT 28 V | 0.25 µm RF HEMT | V_D 28, V_G −2 V | Drain 0…40 V |
| RW: GaN HEMT 650 V (GS66508) | e-mode p-GaN gate: off at V_G = 0, on by 2 V | V_D 1 V | Gate 0…3 V (transfer curve) |
| RW: Thyristor 800 V (BT151) | 4-layer pnpn; the gate sweep shows latching near V_G ≈ 0.65 V | V_AK 2 V, gate 0 | Gate 0…1 V |

For the blocking state of a power switch, set its gate to 0 V and sweep the
drain or collector to the rated voltage. The field distribution and E_max are
meaningful. There is no avalanche model (see Limitations).

**CMOS inverter notes.** The thresholds come out at +0.54 V / −0.54 V
(+0.48 / −0.49 V with the classical models only, against the long-channel
estimate ±0.46 V). Both transistors are one cross-section deep, so they have
equal widths, and V_M sits below VDD/2 (the PMOS is weaker: surface hole
mobility is ~2.5× below the electron one). Real cells size W_p ≈ 2–3 W_n to
centre V_M. A transfer-curve sweep takes 1–3 minutes, because every point
solves the output node iteratively.

**FinFET notes.** The mesh places explicit nodes at the metal surface, at
both faces of the interfacial layer and at the fin surface ("nodes" and
"clear" template keys), so the stack is exactly 0.5 nm SiO2 + 2.5 nm HfO2 on
all three sides. The Structure tab opens on the YZ section through the gate.
For DIBL, set Drain = 0.05 V and sweep again.

**Round and tapered geometry.** The rounded FinFET and the GAA nanowire are
defined with outlines and open on the triangular mesh. Switch them to the
rectangular mesh to see what a staircase does to the same device: with the
default meshes the rounded FinFET's I_on changes by ~25 % and its gate
leakage by ~5×, the nanowire's I_on by ~15 % (the staircase stack is thicker
in some directions and thinner in others). The rounded fin carries ~6 %
less current than the sharp template - a narrower top, and no corner
turning on early.

**Gate-leakage notes.** The two leakage templates are the same transistor
with the same EOT (1.2 nm), so they have the same threshold and inversion
charge; only the dielectric differs. The I-V tab plots the gate current
(linear and semi-log) and lists J_G at ±0.5 and ±1 V. At negative V_G
(accumulation) electrons tunnel from the gate into the substrate and holes
from the accumulation layer into the gate. The gates do not overlap source
and drain: an overlap adds edge tunnelling that dominates at negative V_G.
The high-k result is direct tunnelling only; real stacks also leak through
traps in the HfO2, and measured reductions are 1e2–1e4×.

## How the core works

* **Meshes.** Every equation is assembled with the box method. On the
  rectangular mesh the matrices are 7-point stencils; a prism mesh arrives
  as a general graph (`api3_set_mesh_graph`: symmetric CSR, per edge its
  length and dual face area, per node its volume and position) with the
  material, doping and contact node lists set node by node. The same
  Scharfetter-Gummel, Poisson, MLDA, Lombardi and tunnelling code runs on
  both; the graph path swaps the stencil kernels for CSR ones, the geometric
  multigrid for an aggregation multigrid on the graph, the normal field for a
  least-squares node gradient, the interface walls for ones computed by the
  mesher, and finds tunnelling paths by walking the best-aligned edge.
* **Poisson.** Damped Newton with line search. Each linear step is solved by
  CG with a symmetric aggregation-multigrid preconditioner.
* **Continuity.** Scharfetter–Gummel fluxes with band-modified potentials,
  covering heterojunctions and band-gap narrowing. The system is solved in
  scaled form (u = c/s), so densities spanning 40+ decades keep uniform
  relative accuracy. The solver is BiCGSTAB with a Petrov–Galerkin
  aggregation multigrid (Slotboom-weighted prolongation).
* **Coupling.** The Gummel map runs on the full state (ψ, φₙ, φₚ) with
  Anderson acceleration: depth 12, history stored in float, and depth reduced
  automatically when RAM is short. Carriers too dilute to matter are mixed
  with plain Gummel steps.
* **Maximum principle.** Quasi-Fermi potentials are clamped to the span of
  the terminal voltages. This pins, for example, a floating MOS inversion
  layer exactly.
* **Floating regions.** A floating quasi-neutral region (the n⁻ base of a
  blocking thyristor, an IGBT drift with its gate off) is rescaled after every
  continuity solve so that its exact global current balance holds.
* **Mobility.** Arora low-field mobility, plus Caughey–Thomas velocity
  saturation on every edge. The driving force is the harmonic mean of the
  quasi-Fermi and electrostatic potential drops, which is zero at equilibrium
  and smooth where the two agree.
* **Recombination and band gap.** SRH with Scharfetter τ(N), plus Auger and
  radiative recombination. Klaassen band-gap narrowing.
* **Quantum confinement (MLDA).** At every semiconductor/insulator interface
  (meshed oxide or oxide gate) the density of states of each valley is cut by
  1 − exp(−(z/λ)²), λ = ħ/√(2 m_z kT), averaged over the node's cell; walls
  in several directions multiply (exact for a rectangular corner with
  Boltzmann statistics). Si electrons: Δ2 valleys (m_z 0.916, λ 1.27 nm, 1/3
  of the states) and Δ4 (0.19, 2.79 nm, 2/3); holes: heavy 0.29 (84 %) and
  light 0.20. The correction enters Poisson, the carrier densities and the
  Scharfetter-Gummel band potentials as a per-node ln(N_c,v) shift.
* **Lombardi surface mobility (Si).** 1/μ = 1/μ_bulk + D/μ_ac + D/μ_sr per
  edge, with D = exp(−d/10 nm) (d = distance to the nearest interface) and
  the field normal to the edge from the potential gradient inside the
  semiconductor; Sentaurus parameters. Velocity saturation acts on top.
* **Gate tunnelling.** Paths are found automatically: from every
  semiconductor node next to an insulator, straight through the meshed
  insulator to a gate metal (or through a face gate's layer stack). Along
  each path the oxide potential follows from the series capacitance of the
  layers, and the Tsu-Esaki current ∫T(E)·ln(1+e^((E_F−E)/kT)) dE uses WKB
  transmission through every layer (electrons from the conduction band,
  holes from the valence band, both directions; direct and Fowler-Nordheim
  tunnelling in one formula; over-barrier emission included). The supply of
  an inversion or accumulation layer is fixed by its pressure on the
  interface (the contact theorem kT·n(0) = kT·n(z) + q∫nF dz over a 5 nm
  column), which makes the current independent of the mesh and of MLDA. The
  exchange enters continuity implicitly (a sink linear in the densities of
  the column, the injection as a source), and the gate current closes
  Kirchhoff's law with the other terminals.
* **Continuation.** Adaptive voltage steps with extrapolation. MOS gates are
  ramped before the drain or collector, as in a real measurement. A linear
  solve that still diverges after the ILU-only retry now fails the Gummel
  step (which is then retried shorter) instead of feeding a wrong iterate
  into the loop.
* **Stop and progress.** `api3_request_stop()` (any thread) ends the running
  solve at the next Gummel iteration, with the last converged bias point
  loaded. `api3_get_progress()` reports the ramp fraction, the Gummel
  iteration and the residual.
* **Terminal currents.** Each terminal's current is the Kirchhoff-exact flux
  out of its contact nodes. When the other terminals resolve it better (for
  example a p⁺ anode against a lightly doped cathode), the current is taken
  from the sum of the others instead. Every current carries a round-off
  **resolution floor**; the GUI marks points below it as hollow.

## Limitations

* **No impact ionisation or band-to-band tunnelling.** Avalanche and Zener
  breakdown voltages are not predicted, and there is no GIDL. Blocking sweeps
  give fields and leakage only, and the I-V "knee" marker is a heuristic.
* **No optical generation.** Solar cells and photodiodes are simulated dark.
* **Quantum confinement is MLDA only.** It reproduces the threshold shift and
  the inversion-layer centroid, but there are no subband energies, and at
  very high fields it underestimates the capacitance penalty (0.08 against
  0.1–0.3 nm from Schrödinger-Poisson). Heterostructure 2DEGs (GaN, GaAs)
  stay classical.
* **Gate tunnelling:** no image-force barrier lowering, no trap-assisted
  tunnelling (real high-k stacks leak more than the direct-tunnelling
  value), no valence-electron (EVB) tunnelling from the substrate.
* **No thermionic emission at heterointerfaces.** An abrupt HBT therefore has
  a low β (about 10); real HBTs grade the emitter-base junction.
* **Drift-diffusion transport.** No ballistic or quasi-ballistic transport,
  strain or remote-phonon scattering in high-k stacks: nanoscale on-currents
  (FinFET) are drift-diffusion estimates. Lombardi has Si parameters only;
  SiC, GaN and the other materials keep their bulk mobility at interfaces
  (SiC MOSFET channels are therefore far too good: 20–40 cm²/Vs in reality).
* **Not modelled:** self-heating, traps, incomplete ionisation (Mg in GaN
  enters as its active density), field plates, and transients/AC (the
  inverter transfer curve is DC).
* **The p-GaN gate is ohmic.** A Schottky metal on p-GaN leaves the p-GaN
  floating, with a DC potential set by gate leakage that is not modelled.
* **Low-current accuracy.** Near the resolution floor, currents are good to
  about 0.1–1 %. Well above it, they are good to the tolerances shown above.
* **Speed.** Gummel iteration is slow under strong high-level injection
  (IGBT on-state, latched thyristor) and in saturated channels. Most templates
  solve in seconds; the GaN RF HEMT takes about 1.5 min and the 1200 V IGBT
  about 2–6 min on one core. The Gummel loop stops converging in a few
  high-current regimes:
  * the 1200 V IGBT near V_CE ≈ 2 V;
  * trench MOSFETs at full gate drive beyond V_DS ≈ 3–7 V (5–7 kA/cm²);
  * the superjunction beyond V_DS ≈ 5 V.

  The default sweeps stay inside the converging range. If you push past it,
  the GUI shows the point as not converged and reports the bias it actually
  reached.
* **Memory.** About 760 bytes per node, or about 9 million nodes in 7 GB
  (a prism mesh about 30 % more per node). The GUI refuses meshes larger
  than 80 % of the free RAM.
* **Triangular mesh.** It is a 2-D triangulation extruded along one axis,
  not a general tetrahedral mesh: an outline must be normal to the prism axis
  (others fall back to a staircase), and geometry along the axis stays
  blocky. Outlines are convex (rounded) polygons and circles. Gate tunnelling
  through a curved stack follows each normal path with a planar WKB barrier
  (the field is exact; the film's curvature is ignored in the barrier
  shape). MLDA counts the two nearest walls per node.

## Changes in v7.4

**Triangular (prism) mesh** - a mesh switch on the Device page.

* `semimesh.py`: Delaunay triangulation of the cross-section extruded along
  a prism axis; box geometry keeps the rectangular mesh's nodes near every
  feature and thins them elsewhere; outlines get normal rows with node pairs
  straddling each interface, boundary layers and aligned gate-stack layers.
* **Outlines (shapes)** for regions and box contacts: rounded convex
  polygons and circles, offset families for conformal films, gate electrodes
  wrapped round them. On the rectangular mesh they become staircases.
* **Core:** general box-method graphs (CSR) next to the tensor path, with
  their own aggregation multigrid (strength-based pairwise matching),
  least-squares node gradients, host-supplied interface walls, node-list
  contacts (also usable on the rectangular mesh) and tunnelling paths along
  aligned edges.
* **GUI:** the Structure tab draws outlines and the triangles; the 3D Slicer
  draws prism cross-sections on the triangulation; the Results header shows
  the prism node count; Save/Open keeps outlines and the mesh choice; the
  solution CSV exports the prism nodes.

**Templates:** `FinFET (tapered, rounded fin)` and `GAA nanowire FET
(round)`.

**Validation:** section 13 (prism path = tensor path, thinned-mesh
agreement, round MOS capacitor against radial Poisson-Boltzmann, GAA FET).

## Changes in v7.3

**Physics**

* **Gate tunnelling** (Solver page switch): direct and Fowler-Nordheim
  tunnelling of electrons and holes through any gate stack, meshed (box
  gates) or given as a layer list (face gates). Gate currents appear in the
  Results table, the I-V semi-log panel and the 3D Slicer (log10|J gate|).
* **High-k materials:** HfO2 (k 22, ΔE_c 1.5 eV) and Si3N4 (k 7.5) join SiO2
  and Al2O3; gate stacks on face gates (`SiO2 0.5, HfO2 2.5` in the Contacts
  page) with the EOT computed from the stack.
* **Quantum confinement (MLDA)** at oxide interfaces; the Bands tab draws the
  quantum-corrected band edges Ec + Λn and Ev − Λp.
* **Lombardi surface mobility** for Si inversion and accumulation layers.
* Templates now carry their model switches; NMOS, PMOS, the CMOS inverter
  and the FinFET use all three new models.

**Templates**

* **FinFET** rebuilt with a 0.5 nm SiO2 + 2.5 nm HfO2 stack (EOT 0.94 nm)
  meshed layer by layer; φm 4.50 eV.
* New category "Gate dielectrics & tunnelling": `Gate leakage: 1.2nm SiO2
  NMOS`, `Gate leakage: HfO2 high-k NMOS` (same EOT) and `Fowler-Nordheim:
  8nm SiO2 MOS`.

**Solver and mesh**

* A continuity solve that diverges even after the ILU retry now fails the
  step instead of being used. Before, the 4H-SiC 2 kV validation passed or
  failed depending on round-off (an infinite residual could even be accepted
  as converged); it now passes robustly and faster.
* Mesh: pinned interface nodes replace graded nodes closer than 0.05 nm, and
  templates can clear generated nodes from thin stacks ("clear" key).

**Validation:** sections 11 (tunnelling) and 12 (MLDA, Lombardi) - see above.

## Changes in v7.2

**Interface (rebuilt)**

* **Layout.** Toolbar plus a hideable, resizable editor drawer. The window
  size adapts to the screen (1280×720 and up), and every editor page scrolls
  vertically and horizontally. The mouse wheel works on Windows, macOS and
  Linux/WSLg; Shift + wheel scrolls sideways. Tables keep their own scroll
  bars, so they never widen a page.
* **Text size and settings.** View → Text size (80–175 %) scales the UI and
  the plot labels. Window size, drawer width, text size and the last device
  are remembered.
* **Device library.** Organised in categories, with a browser that shows
  each template's full description.
* **Structure tab.** XY/YZ/XZ cross-sections at any position, drawn with the
  real mesh lines and a contact legend. Hover to read a region's material and
  doping. The 3-D views now show depth downwards.
* **Plots.** Slicer maps are drawn on the true non-uniform mesh (pcolormesh).
  The current map shows the plane containing the cut line.
* **Jobs.** Run and Sweep work in the background with a real progress bar,
  and ■ Stop / Esc stops them. The I-V curve builds up live. When only
  voltages changed, a run continues from the previous solution.
* **Results tab.** A terminal table, extracted parameters and a run log.
* **Help.** A scrollable help window (F1) with topics and search replaces
  the old message box. About is scrollable too.
* **Files.** Save/Open device (JSON) and sweep CSV export.

**Simulation**

* **Nets.** Wired contacts sweep together and have their currents summed.
* **Floating terminals.** A terminal can be solved for zero current at every
  sweep point, which gives transfer curves, with adaptive refinement of steep
  segments.
* **New templates:** CMOS Inverter (0.25 µm) and FinFET (bulk tri-gate).
  NMOS/PMOS now use the common y = 0 top-surface convention.
* **Mesh.** Templates can request explicit nodes and extra boundary layers
  (the FinFET uses them for its 1 nm oxide).
* **Core.** Cooperative stop and a progress record for the GUI; the physics
  is unchanged.
* **Validation.** Section 10 (CMOS thresholds and transfer curve, FinFET SS,
  DIBL and I_on/I_off).

## Changes in v7.1

**Solver**

* Anderson mixing now runs on (ψ, φₙ, φₚ). Mixing ψ alone saw a map with
  hidden state, which stalled saturated channels. HEMT, NMOS and VDMOS solves
  are 2–5× faster, and the GaN HEMT went from never finishing to 18 s.
* Face mobilities are stored in double: float rounding left a 1e-8 V noise
  floor. The velocity-saturation driving force is now smooth.
* New floating-region balance step, and exact quasi-Fermi clamping.
* Noise-floor acceptance.
* Gates are ramped first during continuation.
* Each terminal reports its better-resolved current.
* Built with `-fno-finite-math-only`, so the NaN checks work under
  `-ffast-math`.

**Templates**

* **IGBT:** redesigned. It was pinched off (no JFET implant, 5e13 drift) and
  had no field stop.
* **SiC trench MOSFET:** the shield had wrapped the sidewall and blocked the
  channel exit.
* **n-JFET:** the old 110 nm channel was pinched off at 0 V.
* **GaN e-mode HEMT:** added the missing −σ at the p-GaN/AlGaN interface,
  a realistic barrier, and an ohmic p-GaN gate.
* **HBT:** now vertical, with a 4e19 base; β went from 0.05 to 10.
* **Thyristor:** now demonstrates gate triggering.
* **Power switches:** default to the on-state, not a short-circuit bias.
* **Diode sweeps:** forward is now positive.
* **Mesh:** every oxide face and buried-electrode face now gets a node on
  each side. Before, the far face snapped to the previous grid line: a 100 nm
  trench oxide became about 165 nm on one sidewall, which unbalanced the two
  channels of a trench cell by up to 70 %.
