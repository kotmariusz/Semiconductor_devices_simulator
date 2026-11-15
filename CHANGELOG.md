# Changelog

## 0.2.3 (2026-11-15) - the current animated in 3-D

- Current tab: `▶ Animate` with `3-D` ticked also animates the 3-D view -
  particles travel through the device along the 3-D current (the field of
  the whole device, not one plane), started with probability ∝ |J| and
  moving at a speed that grows with log |J|, measured in the drawn box so a
  thin device does not slow them. The view can be turned while it runs (the
  frame is projected at the present angle over the kept background), and
  untick `3-D arrows` to see the particles alone. `GIF…` saves the 3-D view
  too, as it is turned. The map and the 3-D view can run together or alone.
- Validation: 3-D particles in a helical field keep their distance from the
  axis (5.5e-4 after a quarter turn, in a flattened view box) and move with
  the drift: 240 checks (185 without the templates).

## 0.2.2 (2026-11-15) - current flow pictures, R_DS(on), capacitances

- Current tab: `Show` picks the total current J, the electron flow (−J_n,
  where the electrons go) or the hole flow; `Hide below` drops arrows,
  flowlines and particles where |J| is under a fraction of its maximum
  (leakage paths vanish); `log |J|` colours the map on a log scale.
- Flowlines: lines that each carry the same share of the in-plane current -
  contours of a stream function built from the finite-volume face currents
  of the solution, exact in a 2-D device (path-independent to round-off);
  in a VDMOS cell they show the electrons squeezing through the JFET neck
  and fanning out at about 45° into the drift region. Where current leaves
  the plane (3-D devices) evenly spaced streamlines are drawn and labelled.
- `▶ Animate`: particles flow along the current of the in-plane map, started
  with probability ∝ |J| and moving at a speed that grows with log |J|;
  `Particles`, `Speed`, and `GIF…` saves the animated map.
- R_DS(on): the I-V tab extracts it from an output curve (slope at
  V_DS → 0) or a transfer curve taken at a small V_DS - per simulated cell,
  specific (mΩ·cm², cell footprint) and, with the new `Power device` section
  of the Sweep page (active area, series resistance outside the mesh,
  package), for the die. Analysis → `On-resistance by region…` splits the
  dissipated power J·∇E_F over the regions (channel, JFET, drift,
  substrate) and gives each one's share R_i = P_i/I².
- Capacitances: `Capacitances (quasi-static C-V)` on the Sweep page adds the
  capacitance matrix C_ij = dQ_i/dV_j at every point (each terminal moved by
  a few mV and re-solved). A terminal's charge is Gauss's law around its part
  of the device (its contacts and the semiconductor of the same majority
  carrier connected to them), so no large node charges cancel. A MOSFET
  drain sweep gives the datasheet curves C_iss, C_oss, C_rss (new C-V panel,
  log axis) and Q_oss, E_oss, C_o(tr), C_o(er); a diode or MOS capacitor
  gives C(V). Exact at low frequency in the off state (reciprocity and
  charge conservation are reported, ~1e-6); approximate when a channel
  joins two terminals (flagged). Analysis → `Capacitances at this bias…`
  gives the full matrix. CSV export includes the curves.
- Results: the largest field in the insulators (a gate oxide) and where it
  is - the design limit of SiC trench MOSFETs.
- Defects: `Interface traps 4H-SiC/SiO2: band tail`, acceptor states rising
  exponentially towards E_c (Potbhare et al. 2006/2008: D_edge 5.7e13
  cm⁻²eV⁻¹, E_t 67 meV).
- Limitations documented: the Lombardi surface mobility is silicon only (a
  4H-SiC channel keeps its bulk mobility, so SiC channel resistance comes
  out several times too low), no c-axis mobility anisotropy, uniform
  doping only, capacitances quasi-static (no AC solver, no switching gate
  charge).
- Device files keep the C-V option and the power-device settings (checked
  when a file is opened).
- Validation: section 19 (25 checks: junction and MOS capacitance against
  the depletion formula and the exact low-frequency curve, Gauss oxide
  field, reciprocity and charge conservation of a VDMOS, C_o(tr)/C_o(er),
  R_DS(on) slope and Joule shares of a resistor bar, VDMOS flowlines against
  the drain current, particle tracking, the SiC band tail): 238 checks (183
  without the templates).

## 0.2.1 (2026-11-10) - fixes from the physics and code review

- Material parameters: silicon low-field mobility is now Masetti et al.
  1983 (it mixed Masetti's and Arora's parameters: electrons were 6–14 % too
  high at 1e16–1e17 cm⁻³, 50 % at 1e20); silicon saturation velocity and β
  follow Canali's temperature laws (v_sat,n 22 % lower at 400 K); every
  material has its own mobility temperature exponents (silicon's were used
  for all); GaAs and InP electron mobility from Sotoodeh et al. 2000.
- 4H-SiC: radiative coefficient 1.5e-12 cm³/s (it was 1e5 too small), hole
  mobility N_ref 1.76e19 cm⁻³ (it was 7.5x too small: holes 17–25 % too slow
  at 1e18–1e19), electron affinity 3.7 eV (SiC/SiO₂ offset 2.75 eV).
- Si radiative coefficient 4.73e-15 cm³/s (Trupke 2003); the neutron-damage
  (Perugia) traps take the cross-sections of the cited paper; AlGaN bowing
  cited to the right paper.
- Avalanche: a bias point counts as converged only when the impact-
  ionisation generation has settled (its total and the extrapolated rest of
  the multiplication series) or lies below the current resolution; currents
  just below breakdown crept up by 15–25 % on re-solving, now by a few per
  cent at most, and that remaining uncertainty is added to the resolution
  floor. The breakdown voltages are unchanged.
- I-V analysis: the diode ideality uses kT/q at the sweep's temperature and
  fits the straightest two decades (no series-resistance roll-off); V_th
  subtracts V_DS/2 in the linear region and uses √I_D in saturation; SS is
  the steepest decade above the resolution floor; f_T is taken per interval
  at the right collector current (the Kirk peak sat at the wrong current on
  coarse sweeps); the extraction functions are module-level and tested.
- Results: the convergence flag belongs to the solution shown (pinned and
  A/B results showed the last solve's), and the Kirchhoff line checks the
  contacts' own fluxes (it was zero by construction).
- Robustness: more than 32 contacts or 16 materials is an error (they were
  dropped silently); device files are checked before they are applied (a bad
  file leaves the editor unchanged and names the field); Run refuses
  non-positive sizes, inverted regions and temperatures outside 20–2000 K;
  the core's grid and node getters return NaN for a bad index instead of
  crashing.
- Help and manual: limitations (Boltzmann statistics, full ionisation, no
  negative differential mobility, breakdown above the ionisation data's
  field range, low-current accuracy, limits) and the About text updated.
- Validation: section 18 (25 checks: mobility tables, temperature, limits,
  the extraction functions on exact curves) and a coarse-sweep f_T check:
  213 checks (158 without the templates).

## 0.2 (2026-11-09) - still a test version

- Three new templates:
  - `RW: Si PNP power BJT 250V (MJ21193)` and `RW: Si NPN power BJT 250V
    (MJ21194)`: a 3-D section of a vertical power transistor with a comb
    (interdigitated) emitter - emitter and base-contact fingers joined by
    bus bars; shows current crowding, high injection and quasi-saturation,
    with hFE inside the datasheet limits.
  - `RW: Si lateral PNP LDO 1A (LM2940)`: one round-emitter cell of the
    lateral PNP pass transistor of a bipolar LDO regulator, on the
    triangular mesh; its sweep shows dropout (saturation, rising base and
    substrate current).
- MIT licence ([LICENSE](LICENSE)) and a citation file
  ([CITATION.cff](CITATION.cff)); a note on trademarks and part numbers.
- Sources: the datasheets and patents behind the new templates, Hauser,
  Kirk and Webster for the effects they show, five citations that were
  missing (Ibbetson, Vegard, Fang-Howard, Dowell-Jarratt,
  Sutherland-Hodgman) and updated company names; the Lombardi
  surface-mobility parameters are traced to Lombardi et al. (1988), with an
  open-access paper that lists them.
- Impact ionisation (avalanche breakdown) in the core: Chynoweth law with
  the van Overstraeten - de Man coefficients for Si (GaAs, 4H-SiC and GaN
  from the literature), an optical-phonon temperature factor, F along the
  current; Solver page switch, `log10 G avalanche` in the 3D Slicer, and an
  avalanche breakdown voltage in the I-V tab instead of the old heuristic
  knee (kept for runs with the model off).
- BJT analysis in the I-V tab: every sweep point stores the device's
  electron and hole charge; Gummel sweeps get a third panel with the current
  gain and the quasi-static f_T against I_C (Kirk effect: the fall past the
  f_T peak); output curves get the Early voltage.
- Current tab: `Contact` and the panel 'Current under a contact' - the
  current density through a contact node by node (emitter current
  crowding), with max/mean and centre/edge ratios.
- Comparing designs: `Pin sweep` and `Pins…` (reference curves drawn under
  later sweeps, with their own parameters; show/hide, rename, colour,
  remove, CSV, show the fields of a pin) and `A/B compare…` (capture two
  designs - whole editor state - or two variants of a template, run the
  same sweep on both).
- Visibility: `Panels ▾` on every plot tab hides panels (the rest get the
  room); `Layers…` (Structure and Current tabs) hides regions, contacts and
  defects in the 3-D view (optionally in the section), sets the opacity of
  each part - or of all regions and all contacts (see-through contacts) -
  and can show the insulators in 3-D.
- Current tab: the device structure under the current arrows - `Map`:
  `|J|`, `structure` or `|J| + structure` in the in-plane map, `Structure in
  3-D` around the 3-D arrows, with its own `Opacity` slider (parts hidden or
  given an opacity in `Layers…` keep it).
- I-V tab `Table…`: the parameters of the last sweep and the pinned sweeps
  side by side (V_th, SS, I_on/I_off, n, J_0, V_A, BV, β, f_T, V_M, noise
  margins …) with the temperature and physics models of each, differences
  from a chosen reference column (mV, ×, %), copy as tab-separated text,
  CSV export.
- Five templates with A/B variants (55 in all), group `BJT effects &
  breakdown`: Kirk effect (collector 1e16 / 1e17), emitter current crowding
  (one 20 um emitter / four 5 um fingers), Early effect (base 5e16 / 2e17),
  Si p+n avalanche diode (impact ionisation on / off), NPN BV_CBO / BV_CES /
  BV_CEO.
- A floating terminal whose swept terminal is the one observed (the open
  base of BV_CEO) is plotted as an I-V curve with the open terminal's
  voltage. With impact ionisation on, sweeps end after two failed points in
  a row.
- Device files remember the impact-ionisation switch; the old files open as
  before.
- Validation: section 17 (avalanche breakdown against the ionisation
  integral and Baliga's closed form, its temperature dependence, the stored
  charge, the transit time of a short-base diode, the Early-voltage and f_T
  extraction); section 9 covers the 55 templates (187 checks with
  `--templates`).

## 0.1 (2026-10-08) - first public version (test version)

The first public release. SemiSim 3D was developed privately before this
(internal versions 7.x); 0.1 starts the public numbering.

- 3-D drift-diffusion core in C with OpenMP: box method on graded
  rectangular meshes or triangular prisms, Scharfetter-Gummel fluxes,
  Newton-Poisson with multigrid CG, BiCGSTAB with aggregation multigrid, a
  Gummel map with Anderson acceleration, adaptive bias continuation with a
  sweep predictor, adaptive mesh refinement.
- Physics: heterojunctions, band-gap narrowing, doping- and field-dependent
  mobility with surface mobility, SRH/Auger/radiative recombination, quantum
  confinement (MLDA), gate tunnelling, metal contacts with real barrier
  heights, GaN polarisation, defects (deep levels, interface traps, voids,
  particles, dislocations, radiation damage).
- Python/Tk interface: device editor, 47 templates, bias points and sweeps,
  nets and floating terminals, 3-D structure and solution views, I-V
  analysis, help with the full manual.
- Linux, WSL2 and macOS: `install.sh`, automatic build of the core
  (`semibuild.py`), CI on Ubuntu and macOS.
- Validation: `validate.py`, 167 checks against closed-form theory (47 templates).

Known gaps: see Limitations in [docs/MANUAL.md](docs/MANUAL.md#limitations).
