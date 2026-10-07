# Changelog

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
- Validation: `validate.py`, 167 checks against closed-form theory.

Known gaps: see Limitations in [docs/MANUAL.md](docs/MANUAL.md#limitations).
