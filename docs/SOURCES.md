# Sources

SemiSim 3D is built from published physics. This page lists where its models, parameter values, device structures and test references come from, and the software it uses. The numbers in the program are typical literature values: where sources disagree, one representative value was picked, so a number can differ somewhat from any single reference below. Measured defect properties (capture cross sections above all) scatter by up to 10x or more between studies.

The same list is in the program under Help → Sources.

## Textbooks and data compilations

- S. M. Sze and K. K. Ng, Physics of Semiconductor Devices, 3rd ed., Wiley, 2007 - device physics, Schottky barrier heights, Richardson constants, material properties.
- S. Selberherr, Analysis and Simulation of Semiconductor Devices, Springer, 1984 - drift-diffusion equations, box method, numerical methods.
- Y. Taur and T. H. Ning, Fundamentals of Modern VLSI Devices, 2nd ed., Cambridge University Press, 2009 - MOSFET and bipolar theory used for templates and checks.
- B. J. Baliga, Fundamentals of Power Semiconductor Devices, 2nd ed., Springer, 2019 - power diodes, MOSFETs, IGBTs, JBS diodes, thyristors.
- T. Kimoto and J. A. Cooper, Fundamentals of Silicon Carbide Technology, Wiley-IEEE Press, 2014 - 4H-SiC properties, contacts, defects (Z1/2, EH6/7, stacking faults, micropipes).
- E. H. Nicollian and J. R. Brews, MOS (Metal Oxide Semiconductor) Physics and Technology, Wiley, 1982 - interface traps and oxide charge.
- M. E. Levinshtein, S. L. Rumyantsev and M. S. Shur (eds.), Properties of Advanced Semiconductor Materials: GaN, AlN, InN, BN, SiC, SiGe, Wiley, 2001 - GaN and SiC parameters.
- Ioffe Institute, NSM Archive - Physical Properties of Semiconductors, http://www.ioffe.ru/SVA/NSM/Semicond/ - material parameters of Si, Ge, GaAs, InP, 4H-SiC, GaN and the alloys.
- Synopsys, Sentaurus Device User Guide - the parameter set of the Lombardi surface mobility and the form of the doping-dependent SRH lifetime.

## Transport equations and discretisation

- D. L. Scharfetter and H. K. Gummel, "Large-signal analysis of a silicon Read diode oscillator," IEEE Trans. Electron Devices, vol. 16, no. 1, pp. 64-77, 1969 - the Scharfetter-Gummel current discretisation.
- H. K. Gummel, "A self-consistent iterative scheme for one-dimensional steady state transistor calculations," IEEE Trans. Electron Devices, vol. 11, no. 10, pp. 455-465, 1964 - the Gummel iteration.
- J. W. Slotboom, "Computer-aided two-dimensional analysis of bipolar transistors," IEEE Trans. Electron Devices, vol. 20, no. 8, pp. 669-679, 1973 - Slotboom variables (used to weight the multigrid of the continuity equations).
- R. E. Bank, D. J. Rose and W. Fichtner, "Numerical methods for semiconductor device simulation," IEEE Trans. Electron Devices, vol. 30, no. 9, pp. 1031-1041, 1983 - the box method on triangular (Delaunay/Voronoi) meshes.

## Physical models

### Mobility

- N. D. Arora, J. R. Hauser and D. J. Roulston, "Electron and hole mobilities in silicon as a function of concentration and temperature," IEEE Trans. Electron Devices, vol. 29, no. 2, pp. 292-295, 1982 - doping-dependent mobility (Si parameters).
- D. M. Caughey and R. E. Thomas, "Carrier mobilities in silicon empirically related to doping and field," Proc. IEEE, vol. 55, no. 12, pp. 2192-2193, 1967 - doping and field dependence (velocity saturation).
- C. Canali, G. Majni, R. Minder and G. Ottaviani, "Electron and hole drift velocity measurements in silicon and their empirical relation to electric field and temperature," IEEE Trans. Electron Devices, vol. 22, no. 11, pp. 1045-1047, 1975 - Si saturation velocities and exponents.
- C. Lombardi, S. Manzini, A. Saporito and M. Vanzi, "A physically based mobility model for numerical simulation of nonplanar devices," IEEE Trans. Computer-Aided Design, vol. 7, no. 11, pp. 1164-1171, 1988 - surface mobility.
- M. Roschke and F. Schwierz, "Electron mobility models for 4H, 6H, and 3C SiC," IEEE Trans. Electron Devices, vol. 48, no. 7, pp. 1442-1447, 2001 - 4H-SiC mobility.

### Band structure and statistics

- D. B. M. Klaassen, J. W. Slotboom and H. C. de Graaff, "Unified apparent bandgap narrowing in n- and p-type silicon," Solid-State Electronics, vol. 35, no. 2, pp. 125-129, 1992 - band-gap narrowing.
- C. D. Thurmond, "The standard thermodynamic functions for the formation of electrons and holes in Ge, Si, GaAs, and GaP," J. Electrochem. Soc., vol. 122, no. 8, pp. 1133-1141, 1975 - band gaps of Si, Ge and GaAs versus temperature.
- M. A. Green, "Intrinsic concentration, effective densities of states, and effective mass in silicon," J. Appl. Phys., vol. 67, no. 6, pp. 2944-2954, 1990 - Si effective densities of states.
- I. Vurgaftman, J. R. Meyer and L. R. Ram-Mohan, "Band parameters for III-V compound semiconductors and their alloys," J. Appl. Phys., vol. 89, no. 11, pp. 5815-5875, 2001 - GaAs, InP and GaN band gaps.
- I. Vurgaftman and J. R. Meyer, "Band parameters for nitrogen-containing semiconductors," J. Appl. Phys., vol. 94, no. 6, pp. 3675-3696, 2003 - AlGaN band gap and bowing.
- S. Adachi, "GaAs, AlAs, and AlxGa1-xAs: Material parameters for use in research and device applications," J. Appl. Phys., vol. 58, no. 3, pp. R1-R29, 1985 - AlGaAs band gap and band offsets.
- R. L. Anderson, "Germanium-gallium arsenide heterojunctions," IBM J. Res. Dev., vol. 4, no. 3, pp. 283-287, 1960 - the electron-affinity rule for heterojunctions.
- O. Ambacher et al., "Two-dimensional electron gases induced by spontaneous and piezoelectric polarization charges in N- and Ga-face AlGaN/GaN heterostructures," J. Appl. Phys., vol. 85, no. 6, pp. 3222-3233, 1999 - polarisation charge of the GaN HEMTs.

### Recombination

- W. Shockley and W. T. Read, "Statistics of the recombinations of holes and electrons," Phys. Rev., vol. 87, no. 5, pp. 835-842, 1952, and R. N. Hall, "Electron-hole recombination in germanium," Phys. Rev., vol. 87, no. 2, p. 387, 1952 - SRH recombination and trap occupancy.
- J. G. Fossum and D. S. Lee, "A physical model for the dependence of carrier lifetime on doping density in nondegenerate silicon," Solid-State Electronics, vol. 25, no. 8, pp. 741-747, 1982, and D. J. Roulston, N. D. Arora and S. G. Chamberlain, "Modeling and measurement of minority-carrier lifetime versus doping in diffused layers of n+-p silicon diodes," IEEE Trans. Electron Devices, vol. 29, no. 2, pp. 284-291, 1982 - doping-dependent lifetime.
- J. Dziewior and W. Schmid, "Auger coefficients for highly doped and highly excited silicon," Appl. Phys. Lett., vol. 31, no. 5, pp. 346-348, 1977 - Si Auger coefficients.

### Quantum confinement and tunnelling

- G. Paasch and H. Übensee, "A modified local density approximation. Electron density in inversion layers," Phys. Status Solidi B, vol. 113, no. 1, pp. 165-178, 1982 - the MLDA quantum correction.
- R. Tsu and L. Esaki, "Tunneling in a finite superlattice," Appl. Phys. Lett., vol. 22, no. 11, pp. 562-564, 1973 - the Tsu-Esaki tunnelling current (with WKB transmission).
- R. H. Fowler and L. Nordheim, "Electron emission in intense electric fields," Proc. R. Soc. Lond. A, vol. 119, no. 781, pp. 173-181, 1928.
- J. Robertson, "High dielectric constant oxides," Eur. Phys. J. Appl. Phys., vol. 28, no. 3, pp. 265-291, 2004 - band offsets of HfO2, Al2O3 and Si3N4 on Si.

### Metal contacts

- C. R. Crowell and S. M. Sze, "Current transport in metal-semiconductor barriers," Solid-State Electronics, vol. 9, no. 11-12, pp. 1035-1048, 1966 - thermionic-emission-diffusion and image-force lowering.
- C. R. Crowell, "The Richardson constant for thermionic emission in Schottky barrier diodes," Solid-State Electronics, vol. 8, no. 4, pp. 395-399, 1965, and J. M. Andrews and M. P. Lepselter, "Reverse current-voltage characteristics of metal-silicide Schottky diodes," Solid-State Electronics, vol. 13, no. 7, pp. 1011-1023, 1970 - Richardson constants (n-Si 112, p-Si 32 A/cm²K²).
- F. A. Padovani and R. Stratton, "Field and thermionic-field emission in Schottky barriers," Solid-State Electronics, vol. 9, no. 7, pp. 695-707, 1966.
- C. Y. Chang and S. M. Sze, "Carrier transport across metal-semiconductor barriers," Solid-State Electronics, vol. 13, no. 6, pp. 727-740, 1970, and A. Y. C. Yu, "Electron tunneling and contact resistance of metal-silicon contact barriers," Solid-State Electronics, vol. 13, no. 2, pp. 239-247, 1970 - tunnelling through the barrier, contact resistivity.
- A. M. Cowley and S. M. Sze, "Surface states and barrier height of metal-semiconductor systems," J. Appl. Phys., vol. 36, no. 10, pp. 3212-3220, 1965 - the pinning model for metal/semiconductor pairs without a measured barrier.
- H. B. Michaelson, "The work function of the elements and its periodicity," J. Appl. Phys., vol. 48, no. 11, pp. 4729-4733, 1977 - metal work functions.
- A. C. Schmitz, A. T. Ping, M. A. Khan, Q. Chen, J. W. Yang and I. Adesida, "Metal contacts to n-type GaN," J. Electron. Mater., vol. 27, no. 4, pp. 255-260, 1998 - barrier heights on GaN (on 4H-SiC: Kimoto and Cooper).

## Defects

- A. A. Istratov, H. Hieslmair and E. R. Weber, "Iron and its complexes in silicon," Appl. Phys. A, vol. 69, no. 1, pp. 13-44, 1999 - interstitial iron.
- S. Rein and S. W. Glunz, "Electronic properties of the metastable iron-boron pair in silicon," Appl. Phys. Lett., vol. 82, no. 7, pp. 1054-1056, 2003 - the FeB pair.
- W. M. Bullis, "Properties of gold in silicon," Solid-State Electronics, vol. 9, no. 2, pp. 143-168, 1966 - gold.
- K. Graff, Metal Impurities in Silicon-Device Fabrication, 2nd ed., Springer, 2000 - gold, platinum and other metals in silicon.
- M. Seibt, M. Griess, A. A. Istratov, H. Hedemann, A. Sattler and W. Schröter, "Formation and properties of copper silicide precipitates in silicon," Phys. Status Solidi A, vol. 166, no. 1, pp. 171-182, 1998, and A. A. Istratov and E. R. Weber, "Physics of copper in silicon," J. Electrochem. Soc., vol. 149, no. 1, pp. G21-G30, 2002 - copper and Cu3Si precipitates.
- M. Petasecca, F. Moscatelli, D. Passeri and G. U. Pignatel, "Numerical simulation of radiation damage effects in p-type and n-type FZ silicon detectors," IEEE Trans. Nucl. Sci., vol. 53, no. 5, pp. 2971-2976, 2006 - the "Perugia" neutron-damage trap model.
- T. R. Oldham and F. B. McLean, "Total ionizing dose effects in MOS oxides and devices," IEEE Trans. Nucl. Sci., vol. 50, no. 3, pp. 483-499, 2003 - ionizing-dose damage (oxide charge, interface traps).
- V. V. Afanas'ev, M. Bassler, G. Pensl and M. Schulz, "Intrinsic SiC/SiO2 interface states," Phys. Status Solidi A, vol. 162, no. 1, pp. 321-337, 1997 - SiC/SiO2 interface and near-interface traps.
- J. L. Lyons, A. Janotti and C. G. Van de Walle, "Carbon impurities and the yellow luminescence in GaN," Appl. Phys. Lett., vol. 97, no. 15, 152108, 2010 - the carbon acceptor in GaN.
- A. Y. Polyakov and I.-H. Lee, "Deep traps in GaN-based structures as affecting the performance of GaN devices," Mater. Sci. Eng. R, vol. 94, pp. 1-56, 2015 - Fe and other deep levels, threading dislocations in GaN.
- G. M. Martin, A. Mitonneau and A. Mircea, "Electron traps in bulk and epitaxial GaAs crystals," Electron. Lett., vol. 13, no. 7, pp. 191-193, 1977 - EL2.
- J. Y. W. Seto, "The electrical properties of polycrystalline silicon films," J. Appl. Phys., vol. 46, no. 12, pp. 5247-5254, 1975 - grain boundaries.
- Voids, oxide precipitates, COPs, dislocations in silicon, platinum capture cross sections: representative values from the literature above (no single source); the description of each library defect in the program gives the numbers used.

## Numerical methods

- M. R. Hestenes and E. Stiefel, "Methods of conjugate gradients for solving linear systems," J. Res. Natl. Bur. Stand., vol. 49, no. 6, pp. 409-436, 1952 - conjugate gradients (Poisson).
- H. A. van der Vorst, "Bi-CGSTAB: a fast and smoothly converging variant of Bi-CG for the solution of nonsymmetric linear systems," SIAM J. Sci. Stat. Comput., vol. 13, no. 2, pp. 631-644, 1992 - BiCGSTAB (continuity equations).
- P. Vaněk, J. Mandel and M. Brezina, "Algebraic multigrid by smoothed aggregation for second and fourth order elliptic problems," Computing, vol. 56, no. 3, pp. 179-196, 1996, and Y. Notay, "An aggregation-based algebraic multigrid method," Electron. Trans. Numer. Anal., vol. 37, pp. 123-146, 2010 - the aggregation multigrid preconditioners (pairwise aggregation on the prism mesh).
- D. G. Anderson, "Iterative procedures for nonlinear integral equations," J. ACM, vol. 12, no. 4, pp. 547-560, 1965, and H. F. Walker and P. Ni, "Anderson acceleration for fixed-point iterations," SIAM J. Numer. Anal., vol. 49, no. 4, pp. 1715-1735, 2011 - Anderson acceleration of the Gummel map.
- C. B. Barber, D. P. Dobkin and H. Huhdanpaa, "The Quickhull algorithm for convex hulls," ACM Trans. Math. Softw., vol. 22, no. 4, pp. 469-483, 1996 - Qhull, which computes the Delaunay triangulation of the prism mesh (through matplotlib).

## Validation references (validate.py)

- W. Shockley, "The theory of p-n junctions in semiconductors and p-n junction transistors," Bell Syst. Tech. J., vol. 28, no. 3, pp. 435-489, 1949 - the diode current.
- H. K. Gummel, "Measurement of the number of impurities in the base layer of a transistor," Proc. IRE, vol. 49, no. 4, p. 834, 1961 - the Gummel-number check of the BJT.
- M. Lenzlinger and E. H. Snow, "Fowler-Nordheim tunneling into thermally grown SiO2," J. Appl. Phys., vol. 40, no. 1, pp. 278-283, 1969 - the Fowler-Nordheim slope.
- S.-H. Lo, D. A. Buchanan, Y. Taur and W. Wang, "Quantum-mechanical modeling of electron tunneling current from the inversion layer of ultra-thin-oxide nMOSFET's," IEEE Electron Device Lett., vol. 18, no. 5, pp. 209-211, 1997 - direct tunnelling through thin SiO2.
- S. Takagi, A. Toriumi, M. Iwase and H. Tango, "On the universality of inversion layer mobility in Si MOSFET's: Part I - effects of substrate impurity concentration," IEEE Trans. Electron Devices, vol. 41, no. 12, pp. 2357-2362, 1994 - the universal mobility curve; the fit 540/(1+(E_eff/0.9 MV/cm)^1.85) cm²/Vs used in the check is from K. Chen, H. C. Wann, J. Dunster, P. K. Ko, C. Hu and M. Yoshida, "MOSFET carrier mobility model based on gate oxide thickness, threshold and gate voltages," Solid-State Electronics, vol. 39, no. 10, pp. 1515-1518, 1996.
- F. Stern, "Self-consistent results for n-type Si inversion layers," Phys. Rev. B, vol. 5, no. 12, pp. 4891-4899, 1972; T. Ando, A. B. Fowler and F. Stern, "Electronic properties of two-dimensional systems," Rev. Mod. Phys., vol. 54, no. 2, pp. 437-672, 1982; M. J. van Dort, P. H. Woerlee and A. J. Walker, "A simple model for quantisation effects in heavily-doped silicon MOSFETs at inversion conditions," Solid-State Electronics, vol. 37, no. 3, pp. 411-414, 1994 - inversion-layer centroid and the quantum threshold shift.
- Depletion approximation, MOS capacitor, long-channel MOSFET, CMOS inverter and subthreshold swing: Sze and Ng; Taur and Ning. Crowell-Sze theory for the Schottky diodes (above). The remaining checks compare the core with exact solutions computed by the script itself (radial Poisson-Boltzmann, charge neutrality, direct integration) or with the same device on a finer mesh.

## Device templates

- The templates follow textbook device structures (Sze and Ng; Taur and Ning; Baliga; Kimoto and Cooper) and these papers: D. Hisamoto et al., "FinFET - a self-aligned double-gate MOSFET scalable to 20 nm," IEEE Trans. Electron Devices, vol. 47, no. 12, pp. 2320-2325, 2000 (FinFETs); U. K. Mishra, L. Shen, T. E. Kazior and Y.-F. Wu, "GaN-based RF power devices and amplifiers," Proc. IEEE, vol. 96, no. 2, pp. 287-305, 2008 (GaN HEMTs); T. Fujihira, "Theory of semiconductor superjunction devices," Jpn. J. Appl. Phys., vol. 36, no. 10, pp. 6254-6262, 1997 (superjunction MOSFET); A. W. Blakers, A. Wang, A. M. Milne, J. Zhao and M. A. Green, "22.8% efficient silicon solar cell," Appl. Phys. Lett., vol. 55, no. 13, pp. 1363-1365, 1989 (PERC cell); T. H. Ning, "A perspective on SOI symmetric lateral bipolar transistors for ultra-low-power systems," IEEE J. Electron Devices Soc., vol. 4, no. 5, pp. 227-235, 2016 (SOI lateral NPN); P. R. Gray, P. J. Hurst, S. H. Lewis and R. G. Meyer, Analysis and Design of Analog Integrated Circuits, 5th ed., Wiley, 2009 (lateral PNP of bipolar ICs).
- The "RW:" templates are modelled on the device class and ratings given in the public datasheets of these parts: 1N4148, 1N5819, 1N4733A, BPW34, BC547, IRLZ44N, BT151, Wolfspeed C4D (SiC JBS diodes), Wolfspeed C3M (900 V SiC MOSFETs), Infineon IKW40N120 (1200 V IGBT), Wolfspeed CGH40010 (GaN RF HEMT) and GaN Systems GS66508 (650 V GaN HEMT); the 600 V superjunction MOSFET, the 870 nm AlGaAs IR LED, the 1550 nm Ge photodiode and the PERC cell are generic. Dopings and thicknesses follow published design rules for each voltage class; they are not the real parts' dies, and the simulated currents are per simulated cell.

## Software

- Python (https://www.python.org) with Tkinter and Tcl/Tk (https://www.tcl-lang.org) - the interface.
- C. R. Harris et al., "Array programming with NumPy," Nature, vol. 585, pp. 357-362, 2020 - NumPy (https://numpy.org).
- J. D. Hunter, "Matplotlib: A 2D graphics environment," Comput. Sci. Eng., vol. 9, no. 3, pp. 90-95, 2007 - Matplotlib and mplot3d (https://matplotlib.org).
- OpenMP (https://www.openmp.org), built with GCC (https://gcc.gnu.org) or Clang/LLVM (https://clang.llvm.org) and LLVM's OpenMP runtime libomp; on macOS from Homebrew (https://brew.sh).
- Min Ragan-Kelley, appnope (https://github.com/minrk/appnope) - the way the program turns off macOS App Nap.

## Platform notes

- WSLg (WSL2) sends the clicks of a maximised Linux window to the wrong place: https://github.com/microsoft/wslg/issues/1116 and https://github.com/theitush/giverny/issues/78 - why SemiSim turns maximise into "fill the screen" on WSL.
- Matplotlib's Tk backend needs version 3.10 or newer with Tcl/Tk 9: https://github.com/matplotlib/matplotlib/issues/29126
- Tcl/Tk on macOS (the system Python's Tk 8.5 is too old): https://www.python.org/download/mac/tcltk/

## How it was made

- Much of the code and of this documentation was written with the help of an AI assistant (Claude, by Anthropic). Its physics is checked by validate.py against the theory listed above - not against measurements on real devices - so treat the results with the care described in the README (Status).
