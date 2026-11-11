# Changelog

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
