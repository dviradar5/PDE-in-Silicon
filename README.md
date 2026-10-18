# PDE-in-Silicon

MATLAB research framework for simulating **plasma-dispersion-effect (PDE) focusing in crystalline silicon**, motivated by contactless, high-resolution inspection of integrated circuits. A picosecond Laguerre–Gaussian vortex **pump** creates a ring-shaped distribution of electron–hole pairs. Carrier diffusion and recombination modify the silicon's complex refractive index, forming a transient gradient-index (GRIN)-like lens that focuses a delayed Gaussian **probe**. The main configuration uses **775 nm for both beams** [1, 2].

## Model and workflow

1. **Optical excitation:** `laser.m` and `Helper_Functions/` generate Gaussian and vortex beams; Beer–Lambert attenuation determines absorbed pump energy.
2. **Carrier dynamics:** `Simulation_Functions/FCCDiffusion.m` models photoexcited carrier generation, diffusion and recombination (assuming one electron–hole pair per absorbed pump photon).
3. **PDE:** `Simulation_Functions/complexRefractiveIndex.m` and `Helper_Functions/Drude.m` convert carrier concentration into refractive-index and absorption changes.
4. **Probe propagation:** `Simulation_Functions/propagationBPM_rz.m` uses an axisymmetric, scalar, paraxial finite-difference beam-propagation method (FD-BPM) with Crank–Nicolson stepping.
5. **Analysis:** `Tester_Functions/` evaluates probe full width at half maximum (FWHM), focusing and optical losses; `Plot_Functions/` generates figures.

## Run

Key files: `systemParameters.m` (material and numerical settings), `laser.m` (beam definitions), `mainProjectSimulation.m` (main driver), `asaf.m` (experimental comparison), and `Figures_and_Results/` (saved results). `Papers/` contains background literature; `unused/` contains historical implementations.

## Main findings

The final project [1] reports that the simulation reproduces major experimental focusing trends [2, 3]: higher pump energy or a narrower vortex pump generally reduces probe FWHM, while **increased free-carrier absorption reduces transmitted power**. Strong focusing can move the minimum-width plane inside the sample; very narrow pumps can produce multiple axial focusing regions. FWHM reduction alone does **not** establish imaging resolution or contrast: those require target-based evaluation. In the 2024 experiment [2], a 35 nJ pump reduced the measured probe FWHM from approximately **2.8 to 1.2 µm**.

## Scope and limitations

The model assumes radial symmetry, scalar paraxial propagation, idealized beams, simplified carrier generation and free-carrier optical response. It does not fully account for **Auger recombination, band filling, nonlinear absorption, thermal effects, surface imperfections, spherical aberration, multiple interface reflections, or experimental noise**. At high carrier densities or pump energies, these omissions and silicon damage thresholds can limit quantitative predictions [1, 2, 4]. Check grid/domain convergence, pump-energy normalization, transmission and the definition of pump diameter before comparing against experiments. The model is intended for research and parameter exploration, not validated device design.

## References

1. **D. Adar**, *Implementation of an Innovative Photonic Super-Resolution Method for Fault Monitoring in Silicon Microelectronic Chips* (original Hebrew title: *מימוש שיטה פוטונית חדשנית של סופר רזולוציה לניטור תקלות בשבבי סיליקון במיקרו אלקטרוניקה*), B.Sc. final project report, Faculty of Engineering, Bar-Ilan University, October 2026. Supervisors: M. Sinvani and Z. Zalevsky. [Final project report](https://github.com/dviradar5/PDE-in-Silicon/blob/main/Papers/%D7%93%D7%95%D7%97%20%D7%9E%D7%A1%D7%9B%D7%9D.pdf)
2. N. Shabairou, Z. Zalevsky and M. Sinvani, “Laser Beam Self-Focusing in Silicon at an Absorbed Wavelength by a Vortex Beam in the Same Wavelength,” *ACS Omega* **9** (2024), 16969–16975. https://doi.org/10.1021/acsomega.3c08325
3. N. Shabairou, M. Tiferet, Z. Zalevsky and M. Sinvani, “Dynamics of Laser-Induced Tunable Focusing in Silicon,” *Scientific Reports* **12** (2022), 6342. https://doi.org/10.1038/s41598-022-10112-3
4. R. A. Soref and B. R. Bennett, “Electrooptical Effects in Silicon,” *IEEE Journal of Quantum Electronics* **23** (1987), 123–129.
5. H. Pinhas *et al.*, “Plasma Dispersion Effect Based Super-Resolved Imaging in Silicon,” *Optics Express* **26** (2018), 25370–25380. https://doi.org/10.1364/OE.26.025370
