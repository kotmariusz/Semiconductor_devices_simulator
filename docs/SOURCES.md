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
- Synopsys, Sentaurus Device User Guide - the T/300 K form in which the Lombardi surface-mobility parameters are written (the values are those of Lombardi et al. 1988, under Mobility) and the Scharfetter relation for the doping-dependent SRH lifetime (N_ref = 1e16 cm⁻³; the form is that of Fossum and Lee 1982, under Recombination).

## Transport equations and discretisation

- D. L. Scharfetter and H. K. Gummel, "Large-signal analysis of a silicon Read diode oscillator," IEEE Trans. Electron Devices, vol. 16, no. 1, pp. 64-77, 1969 - the Scharfetter-Gummel current discretisation.
- H. K. Gummel, "A self-consistent iterative scheme for one-dimensional steady state transistor calculations," IEEE Trans. Electron Devices, vol. 11, no. 10, pp. 455-465, 1964 - the Gummel iteration.
- J. W. Slotboom, "Computer-aided two-dimensional analysis of bipolar transistors," IEEE Trans. Electron Devices, vol. 20, no. 8, pp. 669-679, 1973 - Slotboom variables (used to weight the multigrid of the continuity equations).
- R. E. Bank, D. J. Rose and W. Fichtner, "Numerical methods for semiconductor device simulation," IEEE Trans. Electron Devices, vol. 30, no. 9, pp. 1031-1041, 1983 - the box method on triangular (Delaunay/Voronoi) meshes.

## Physical models

### Mobility

- G. Masetti, M. Severi and S. Solmi, "Modeling of carrier mobility against carrier concentration in arsenic-, phosphorus-, and boron-doped silicon," IEEE Trans. Electron Devices, vol. 30, no. 7, pp. 764-769, 1983 - Si low-field mobility (arsenic set for electrons, boron set for holes, as the Sentaurus defaults; lattice term scaled by (T/300)^-2.5 / ^-2.2).
- N. D. Arora, J. R. Hauser and D. J. Roulston, "Electron and hole mobilities in silicon as a function of concentration and temperature," IEEE Trans. Electron Devices, vol. 29, no. 2, pp. 292-295, 1982 - the doping- and temperature-dependent (Caughey-Thomas) form used for the other materials; until v0.2 the Si mobility mixed this paper's and Masetti's parameters.
- D. M. Caughey and R. E. Thomas, "Carrier mobilities in silicon empirically related to doping and field," Proc. IEEE, vol. 55, no. 12, pp. 2192-2193, 1967 - doping and field dependence (velocity saturation).
- C. Canali, G. Majni, R. Minder and G. Ottaviani, "Electron and hole drift velocity measurements in silicon and their empirical relation to electric field and temperature," IEEE Trans. Electron Devices, vol. 22, no. 11, pp. 1045-1047, 1975 - Si saturation velocities and exponents, with their temperature laws (v_sat ∝ T^-0.87 / T^-0.52, β ∝ T^0.66 / T^0.17).
- C. Lombardi, S. Manzini, A. Saporito and M. Vanzi, "A physically based mobility model for numerical simulation of nonplanar devices," IEEE Trans. Computer-Aided Design, vol. 7, no. 11, pp. 1164-1171, 1988, doi:10.1109/43.9186 - surface mobility. SemiSim uses the 1988 parameters for Si with the temperature written as T/300 K, so C = C(1988)/300 K (1.74e5/300 = 580 for electrons, 8.84e5/300 = 2947 for holes). The electron values are also listed in the open-access paper A. Ciprut, A. Chelly and A. Karsenty, "Usage and limitation of standard mobility models for TCAD simulation of nanoscaled FD-SOI MOSFETs," Active and Passive Electronic Components, vol. 2015, Art. 460416, 2015, doi:10.1155/2015/460416.
- M. Roschke and F. Schwierz, "Electron mobility models for 4H, 6H, and 3C SiC," IEEE Trans. Electron Devices, vol. 48, no. 7, pp. 1442-1447, 2001 - 4H-SiC electron mobility and its temperature exponents.
- W. J. Schaffer, G. H. Negley, K. G. Irvine and J. W. Palmour, "Conductivity anisotropy in epitaxial 6H and 4H SiC," MRS Proc., vol. 339, p. 595, 1994 - 4H-SiC hole mobility (124, 15.9 cm²/Vs, N_ref 1.76e19 cm⁻³, α 0.34), as tabulated by T. Ayalew, PhD thesis, TU Wien, 2004.
- M. Sotoodeh, A. H. Khalid and A. A. Rezazadeh, "Empirical low-field mobility model for III-V compounds applicable in device simulation codes," J. Appl. Phys., vol. 87, no. 6, pp. 2890-2900, 2000 - GaAs and InP electron mobility and their temperature exponents.
- S. Farahmand et al., "Monte Carlo simulation of electron transport in the III-nitride wurtzite phase materials system: binaries and ternaries," IEEE Trans. Electron Devices, vol. 48, no. 3, pp. 535-542, 2001 - temperature exponents of the GaN mobility.

### Band structure and statistics

- D. B. M. Klaassen, J. W. Slotboom and H. C. de Graaff, "Unified apparent bandgap narrowing in n- and p-type silicon," Solid-State Electronics, vol. 35, no. 2, pp. 125-129, 1992 - band-gap narrowing.
- C. D. Thurmond, "The standard thermodynamic functions for the formation of electrons and holes in Ge, Si, GaAs, and GaP," J. Electrochem. Soc., vol. 122, no. 8, pp. 1133-1141, 1975 - band gaps of Si, Ge and GaAs versus temperature.
- M. A. Green, "Intrinsic concentration, effective densities of states, and effective mass in silicon," J. Appl. Phys., vol. 67, no. 6, pp. 2944-2954, 1990 - Si effective densities of states.
- I. Vurgaftman, J. R. Meyer and L. R. Ram-Mohan, "Band parameters for III-V compound semiconductors and their alloys," J. Appl. Phys., vol. 89, no. 11, pp. 5815-5875, 2001 - GaAs, InP and GaN band gaps.
- I. Vurgaftman and J. R. Meyer, "Band parameters for nitrogen-containing semiconductors," J. Appl. Phys., vol. 94, no. 6, pp. 3675-3696, 2003 - the nitride update; it recommends an AlGaN bowing of 0.7 eV, while the program keeps the 1.0 eV of the 2001 paper (above).
- V. V. Afanas'ev, M. Bassler, G. Pensl, M. J. Schulz and E. Stein von Kamienski, "Band offsets and electronic structure of SiC/SiO2 interfaces," J. Appl. Phys., vol. 79, no. 6, pp. 3108-3114, 1996 - the SiC/SiO2 band offsets that set the 4H-SiC electron affinity (3.7 eV: conduction-band offset 2.75 eV).
- S. Adachi, "GaAs, AlAs, and AlxGa1-xAs: Material parameters for use in research and device applications," J. Appl. Phys., vol. 58, no. 3, pp. R1-R29, 1985 - AlGaAs band gap and band offsets.
- R. L. Anderson, "Germanium-gallium arsenide heterojunctions," IBM J. Res. Dev., vol. 4, no. 3, pp. 283-287, 1960 - the electron-affinity rule for heterojunctions.
- O. Ambacher et al., "Two-dimensional electron gases induced by spontaneous and piezoelectric polarization charges in N- and Ga-face AlGaN/GaN heterostructures," J. Appl. Phys., vol. 85, no. 6, pp. 3222-3233, 1999 - polarisation charge of the GaN HEMTs.
- J. P. Ibbetson, P. T. Fini, K. D. Ness, S. P. DenBaars, J. S. Speck and U. K. Mishra, "Polarization effects, surface states, and the source of electrons in AlGaN/GaN heterostructure field effect transistors," Appl. Phys. Lett., vol. 77, no. 2, pp. 250-252, 2000 - surface donors as the source of the channel electrons (GaN HEMT templates).
- L. Vegard, "Die Konstitution der Mischkristalle und die Raumfüllung der Atome," Z. Phys., vol. 5, no. 1, pp. 17-26, 1921 - linear interpolation of alloy parameters (with the bowing terms of the papers above).

### Recombination

- W. Shockley and W. T. Read, "Statistics of the recombinations of holes and electrons," Phys. Rev., vol. 87, no. 5, pp. 835-842, 1952, and R. N. Hall, "Electron-hole recombination in germanium," Phys. Rev., vol. 87, no. 2, p. 387, 1952 - SRH recombination and trap occupancy.
- J. G. Fossum and D. S. Lee, "A physical model for the dependence of carrier lifetime on doping density in nondegenerate silicon," Solid-State Electronics, vol. 25, no. 8, pp. 741-747, 1982, and D. J. Roulston, N. D. Arora and S. G. Chamberlain, "Modeling and measurement of minority-carrier lifetime versus doping in diffused layers of n+-p silicon diodes," IEEE Trans. Electron Devices, vol. 29, no. 2, pp. 284-291, 1982 - doping-dependent lifetime.
- J. Dziewior and W. Schmid, "Auger coefficients for highly doped and highly excited silicon," Appl. Phys. Lett., vol. 31, no. 5, pp. 346-348, 1977 - Si Auger coefficients.
- T. Trupke et al., "Temperature dependence of the radiative recombination coefficient of intrinsic crystalline silicon," J. Appl. Phys., vol. 94, no. 8, pp. 4930-4937, 2003 - Si radiative coefficient 4.73e-15 cm³/s.
- A. Galeckas, J. Linnros, V. Grivickas, U. Lindefelt and C. Hallin, "Auger recombination in 4H-SiC: unusual temperature behavior," Appl. Phys. Lett., vol. 71, no. 22, pp. 3269-3271, 1997 - 4H-SiC Auger coefficients and radiative coefficient (1.5e-12 cm³/s, as listed by the Ioffe NSM archive).

### Impact ionisation (avalanche breakdown)

- A. G. Chynoweth, "Ionization rates for electrons and holes in silicon," Phys. Rev., vol. 109, no. 5, pp. 1537-1540, 1958 - the field dependence α = a exp(−b/F).
- R. van Overstraeten and H. de Man, "Measurement of the ionization rates in diffused silicon p-n junctions," Solid-State Electronics, vol. 13, no. 5, pp. 583-608, 1970 - the Si coefficients (holes in two field ranges).
- C. R. Crowell and S. M. Sze, "Temperature dependence of avalanche multiplication in semiconductors," Appl. Phys. Lett., vol. 9, no. 6, pp. 242-244, 1966 - the temperature dependence through the optical-phonon mean free path, λ ∝ tanh(ħω/2kT) (Si: ħω = 63 meV), behind the factor γ = tanh(ħω/2kT₀)/tanh(ħω/2kT).
- G. E. Bulman, V. M. Robbins, K. F. Brennan, K. Hess and G. E. Stillman, "Experimental determination of impact ionization coefficients in (100) GaAs," IEEE Electron Device Lett., vol. 4, no. 6, pp. 181-185, 1983 - GaAs; the fit α = a exp(−(b/F)^m) used here as tabulated by R. Quay, Analysis and Simulation of High Electron Mobility Transistors, PhD thesis, TU Wien, 2001 (Table 3.26).
- T. Hatakeyama, T. Watanabe, T. Shinohe, K. Kojima, K. Arai and N. Sano, "Impact ionization coefficients of 4H silicon carbide," Appl. Phys. Lett., vol. 85, no. 8, pp. 1380-1382, 2004 - 4H-SiC (field along the c axis).
- D. Ji, B. Ercan and S. Chowdhury, "Experimental determination of impact ionization coefficients of electrons and holes in gallium nitride using homojunction structures," Appl. Phys. Lett., vol. 115, no. 7, 073503, 2019 - GaN.
- S. L. Miller, "Avalanche breakdown in germanium," Phys. Rev., vol. 99, no. 4, pp. 1234-1241, 1955 - the multiplication law M = 1/(1 − (V/BV)^n), and with it BV_CEO ≈ BV_CBO/β^(1/n) (the open-base template).

### Bipolar transistor analysis (I-V tab)

- A. Ortiz-Conde, F. J. García Sánchez, J. J. Liou, A. Cerdeira, M. Estrada and Y. Yue, "A review of recent MOSFET threshold voltage extraction methods," Microelectron. Reliab., vol. 42, no. 4-5, pp. 583-596, 2002 - V_th by linear extrapolation at maximum g_m (minus V_DS/2) and by the √I_D extrapolation in saturation.

- J. M. Early, "Effects of space-charge layer widening in junction transistors," Proc. IRE, vol. 40, no. 11, pp. 1401-1406, 1952 - base-width modulation and the Early voltage.
- The quasi-static f_T = g_m / (2π dQ/dV_BE) is the charge-control description of the transistor (Taur and Ning, chapters on bipolar transistor performance); Kirk (below, Device templates) for the f_T fall-off at high current, Hauser for emitter current crowding.

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
- D. Macdonald, A. Cuevas and J. Wong-Leung, "Capture cross sections of the acceptor level of iron-boron pairs in p-type silicon by injection-level dependent lifetime measurements," J. Appl. Phys., vol. 89, no. 12, pp. 7932-7939, 2001, doi:10.1063/1.1372156, and S. Rein and S. W. Glunz, "Electronic properties of interstitial iron and iron-boron pairs determined by means of advanced lifetime spectroscopy," J. Appl. Phys., vol. 98, no. 11, 113711, 2005, doi:10.1063/1.2106017 - the FeB pair. The level is E_C - 0.26 ± 0.03 eV in Istratov et al. (above); SemiSim uses E_C - 0.23 eV, as Macdonald et al. do, and representative capture cross sections (σn 5e-15, σp 3e-15 cm²) - measured values scatter widely between these papers.
- W. M. Bullis, "Properties of gold in silicon," Solid-State Electronics, vol. 9, no. 2, pp. 143-168, 1966 - gold.
- K. Graff, Metal Impurities in Silicon-Device Fabrication, 2nd ed., Springer, 2000 - gold, platinum and other metals in silicon.
- M. Seibt, M. Griess, A. A. Istratov, H. Hedemann, A. Sattler and W. Schröter, "Formation and properties of copper silicide precipitates in silicon," Phys. Status Solidi A, vol. 166, no. 1, pp. 171-182, 1998, and A. A. Istratov and E. R. Weber, "Physics of copper in silicon," J. Electrochem. Soc., vol. 149, no. 1, pp. G21-G30, 2002 - copper and Cu3Si precipitates.
- M. Petasecca, F. Moscatelli, D. Passeri and G. U. Pignatel, "Numerical simulation of radiation damage effects in p-type and n-type FZ silicon detectors," IEEE Trans. Nucl. Sci., vol. 53, no. 5, pp. 2971-2976, 2006 - the "Perugia" neutron-damage trap model (levels, cross-sections and introduction rates of its p-type set, as tabulated by F. Moscatelli et al., arXiv:1611.10132, which also lists a modified set).
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
- M. Dowell and P. Jarratt, "A modified regula falsi method for computing the root of an equation," BIT, vol. 11, no. 2, pp. 168-174, 1971 - the Illinois method, which finds the voltage of floating terminals.
- I. E. Sutherland and G. W. Hodgman, "Reentrant polygon clipping," Commun. ACM, vol. 17, no. 1, pp. 32-42, 1974 - clipping of round and sloped outlines in the structure views.

## Validation references (validate.py)

- W. Shockley, "The theory of p-n junctions in semiconductors and p-n junction transistors," Bell Syst. Tech. J., vol. 28, no. 3, pp. 435-489, 1949 - the diode current.
- H. K. Gummel, "Measurement of the number of impurities in the base layer of a transistor," Proc. IRE, vol. 49, no. 4, p. 834, 1961 - the Gummel-number check of the BJT.
- M. Lenzlinger and E. H. Snow, "Fowler-Nordheim tunneling into thermally grown SiO2," J. Appl. Phys., vol. 40, no. 1, pp. 278-283, 1969 - the Fowler-Nordheim slope.
- S.-H. Lo, D. A. Buchanan, Y. Taur and W. Wang, "Quantum-mechanical modeling of electron tunneling current from the inversion layer of ultra-thin-oxide nMOSFET's," IEEE Electron Device Lett., vol. 18, no. 5, pp. 209-211, 1997 - direct tunnelling through thin SiO2.
- S. Takagi, A. Toriumi, M. Iwase and H. Tango, "On the universality of inversion layer mobility in Si MOSFET's: Part I - effects of substrate impurity concentration," IEEE Trans. Electron Devices, vol. 41, no. 12, pp. 2357-2362, 1994 - the universal mobility curve; the fit 540/(1+(E_eff/0.9 MV/cm)^1.85) cm²/Vs used in the check is from K. Chen, H. C. Wann, J. Dunster, P. K. Ko, C. Hu and M. Yoshida, "MOSFET carrier mobility model based on gate oxide thickness, threshold and gate voltages," Solid-State Electronics, vol. 39, no. 10, pp. 1515-1518, 1996.
- F. Stern, "Self-consistent results for n-type Si inversion layers," Phys. Rev. B, vol. 5, no. 12, pp. 4891-4899, 1972; T. Ando, A. B. Fowler and F. Stern, "Electronic properties of two-dimensional systems," Rev. Mod. Phys., vol. 54, no. 2, pp. 437-672, 1982; M. J. van Dort, P. H. Woerlee and A. J. Walker, "A simple model for quantisation effects in heavily-doped silicon MOSFETs at inversion conditions," Solid-State Electronics, vol. 37, no. 3, pp. 411-414, 1994 - inversion-layer centroid and the quantum threshold shift; F. F. Fang and W. E. Howard, "Negative field-effect mobility on (100) Si surfaces," Phys. Rev. Lett., vol. 16, no. 18, pp. 797-799, 1966 - the Fang-Howard estimate of the centroid.
- Avalanche breakdown (section 17): the ionisation-integral criterion and Baliga's closed form BV = 5.34e13 N^(−3/4) V for abrupt Si junctions (Sze and Ng; Baliga); the stored charge of a short-base diode, Q = I W²/(2D) (Sze and Ng).
- Depletion approximation, MOS capacitor, long-channel MOSFET, CMOS inverter and subthreshold swing: Sze and Ng; Taur and Ning. Crowell-Sze theory for the Schottky diodes (above). The remaining checks compare the core with exact solutions computed by the script itself (radial Poisson-Boltzmann, charge neutrality, direct integration) or with the same device on a finer mesh.

## Device templates

- The templates follow textbook device structures (Sze and Ng; Taur and Ning; Baliga; Kimoto and Cooper) and these papers: D. Hisamoto et al., "FinFET - a self-aligned double-gate MOSFET scalable to 20 nm," IEEE Trans. Electron Devices, vol. 47, no. 12, pp. 2320-2325, 2000 (FinFETs); U. K. Mishra, L. Shen, T. E. Kazior and Y.-F. Wu, "GaN-based RF power devices and amplifiers," Proc. IEEE, vol. 96, no. 2, pp. 287-305, 2008 (GaN HEMTs); T. Fujihira, "Theory of semiconductor superjunction devices," Jpn. J. Appl. Phys., vol. 36, no. 10, pp. 6254-6262, 1997 (superjunction MOSFET); A. W. Blakers, A. Wang, A. M. Milne, J. Zhao and M. A. Green, "22.8% efficient silicon solar cell," Appl. Phys. Lett., vol. 55, no. 13, pp. 1363-1365, 1989 (PERC cell); T. H. Ning, "A perspective on SOI symmetric lateral bipolar transistors for ultra-low-power systems," IEEE J. Electron Devices Soc., vol. 4, no. 5, pp. 227-235, 2016 (SOI lateral NPN); P. R. Gray, P. J. Hurst, S. H. Lewis and R. G. Meyer, Analysis and Design of Analog Integrated Circuits, 5th ed., Wiley, 2009 (lateral PNP of bipolar ICs).
- J. R. Hauser, "The effects of distributed base potential on emitter-current injection density and effective base resistance for stripe transistor geometries," IEEE Trans. Electron Devices, vol. 11, no. 5, pp. 238-242, 1964 - emitter current crowding, the reason power transistors have comb (interdigitated) emitters; C. T. Kirk, "A theory of transistor cutoff frequency (fT) falloff at high current densities," IRE Trans. Electron Devices, vol. 9, no. 2, pp. 164-174, 1962 - base push-out (Kirk effect) and quasi-saturation in the lightly doped collector; W. M. Webster, "On the variation of junction-transistor current-amplification factor with emitter current," Proc. IRE, vol. 42, no. 6, pp. 914-920, 1954 - the fall of the current gain at high injection (power BJT and lateral PNP templates).
- J. R. Murkland and J. S. Congdon, "Lateral PNP power transistor," U.S. Patent 4,417,265, 1983, and T. J. Skovmand, "Quiescent current reduction in low dropout voltage regulators," U.S. Patent 4,613,809, 1986 (both National Semiconductor) - the many-emitter lateral PNP pass transistor of bipolar LDO regulators: 0.5 mil round emitters, each ringed by a collector, over an n+ buried layer; about 3 mA per emitter, 350 emitters for 1 A at 0.5 V dropout (the LDO lateral PNP template). National Semiconductor, LM2940/LM2940C 1A Low Dropout Regulator data sheet, DS008822, 2000 - dropout (0.5 V at 1 A) and ground current (30 mA typical at 1 A) the template is compared with; C. Simpson, "Linear and switching voltage regulator fundamentals, Part 1," National Semiconductor (Texas Instruments SNVA558) - how dropout and ground current of an LDO depend on its PNP pass transistor.
- The "RW:" templates are modelled on the device class and ratings given in the public datasheets of these parts: 1N4148, 1N5819, 1N4733A, BPW34, BC547, MJ21193 and MJ21194 (onsemi 250 V / 16 A audio power BJTs; onsemi builds them with a "perforated emitter", the templates with the classic comb), LM2940 (1 A bipolar LDO regulator, National Semiconductor, now Texas Instruments), IRLZ44N, BT151, Wolfspeed C4D (SiC JBS diodes), Wolfspeed C3M (900 V SiC MOSFETs), Infineon IKW40N120 (1200 V IGBT), CGH40010 (GaN RF HEMT; Cree/Wolfspeed, sold by MACOM since 2023) and GaN Systems GS66508 (650 V GaN HEMT; GaN Systems has been part of Infineon since 2023); the 600 V superjunction MOSFET, the 870 nm AlGaAs IR LED, the 1550 nm Ge photodiode and the PERC cell are generic. Dopings and thicknesses follow published design rules for each voltage class; they are not the real parts' dies, and the simulated currents are per simulated cell.
- Part numbers, company and product names (including Sentaurus and Synopsys above) are used only to identify the parts and documents meant; they are trademarks of their owners, who are not connected with SemiSim 3D.

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
