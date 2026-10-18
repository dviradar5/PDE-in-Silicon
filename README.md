# PDE-in-Silicon

MATLAB numerical framework for studying **plasma-dispersion-induced focusing in crystalline silicon**. A picosecond Laguerre–Gaussian vortex pump generates electron–hole pairs, which diffuse and recombine, producing a transient complex-index distribution. A delayed Gaussian probe is propagated through that distribution using an axisymmetric finite-difference beam-propagation method (FD-BPM).

The code supports studies of pump energy, beam waist and waist position, temporal delay, free-carrier concentration, probe propagation, absorption, and exit-plane full width at half maximum (FWHM). This is a research simulation, not a validated manufacturing design tool.

## Physical workflow

1. `laser.m` constructs pump and probe fields with `Helper_Functions/LG0B.m` and `GB.m`.
2. The pump intensity is attenuated through silicon using Beer–Lambert absorption (`laser.intensityProfileBLDumped`).
3. `Simulation_Functions/FCCDiffusion.m` computes the photoinduced carrier density, incorporating generation, cylindrical diffusion and recombination.
4. `Simulation_Functions/complexRefractiveIndex.m` and `Helper_Functions/Drude.m` map the carrier density onto the complex optical index.
5. `Simulation_Functions/propagationBPM_rz.m` propagates the probe through the resulting index map using a radial cylindrical Laplacian, Crank–Nicolson diffraction steps and a radial absorbing mask.
6. `Tester_Functions/FWHM.m` estimates the central-lobe FWHM of the probe.

The model includes diffraction, linear pump absorption, carrier diffusion and a free-carrier Drude contribution. It does **not** fully model high-density band filling, Auger recombination, experimental wavefront aberrations, vector polarization, internal multiple reflections, nonlinear absorption or material damage.

## Repository layout

| Location | Purpose |
|---|---|
| `mainProjectSimulation.m` | Main interactive example / original project driver |
| `asaf.m` | Configuration based on Asaf's pump-compression experiment |
| `systemParameters.m` | Constants, silicon optical properties and model defaults |
| `laser.m` | Gaussian and vortex optical fields and pulse energy |
| `Simulation_Functions/` | FCC generation/diffusion, PDE index and probe propagation |
| `Helper_Functions/` | Beam profiles, Drude model, diameter and plotting conversions |
| `Tester_Functions/` | FWHM, focusing checks, power loss and threshold helpers |
| `Plot_Functions/` | Research plots and visualizations |
| `Figures_and_Results/` | Previously saved experimental/numerical figures |
| `Papers/` | Source literature |
| `unused/` | Historical and experimental implementations |

## Requirements

- MATLAB with the functions and toolboxes used in the chosen script (test on your local MATLAB installation; no specific MATLAB release has been verified here).
- No Python environment is required for the MATLAB simulation.
- MATLAB path includes the repository and helper folders.

```matlab
cd('path/to/PDE-in-Silicon')
addpath(genpath(pwd))
mainProjectSimulation
% Or run: asaf
```

## Pump-width experiment: Figure 7 extension

Copy `pump_width_sweep.m`, `pump_width_sweep_config.m` and `plot_pump_width_sweep.m` into this repository's root, then run:

```matlab
addpath(genpath(pwd));
cfg = pump_width_sweep_config();
results = pump_width_sweep(cfg);
```

The extension compares **20 nJ (Asaf)** and **35 nJ (Nadav, 2024)**, each with four incident probe Gaussian 1/e² waist radii: 1.875, 5, 10 and 100 µm. The last is a **finite-width approximation to a plane wave**, not an exact plane wave. The driver calculates the induced medium once per pump diameter and reuses it for each probe case. Results are saved to `.mat`, `.csv`, `.fig` and `.png` files in `pump_width_results/`.

### Definition of pump width

The horizontal axis is the pump's **peak-to-peak ring diameter on the back surface**. For LG₀¹ intensity,

\[
I(r,z)\propto \left(\frac{2r^2}{w(z)^2}\right)e^{-2r^2/w(z)^2},\qquad
r_{\mathrm{peak}}(z)=\frac{w(z)}{\sqrt{2}},\quad
D_{\mathrm{peak}}(z)=\sqrt{2}\,w(z).
\]

This is **not** the same as the dark central-hole width, which is reported as 0.91 µm in Asaf's measured pump profile. The Figure 7 sweep uses peak-to-peak diameter; convert the measured pump profile using the same metric before putting Asaf's experimental **X** on that axis. The driver leaves the unknown coordinate unset rather than making up a value. Nadav's 2024 article reports an approximately 6 µm ring and a probe FWHM of approximately 1.2 µm (35 nJ), so that reference can be marked directly, subject to confirmation of its diameter metric.

### Reference cases

| Parameter | Asaf experiment | Nadav (2024) |
|---|---:|---:|
| Pump and probe wavelength | 775 nm | 775 nm |
| Pump pulse FWHM | 30 ps | 30 ps |
| Pump energy before interface | 20 nJ | 35 nJ |
| Probe delay | ≈100 ps | 60 ps |
| Sample thickness | ≈20 µm | ≈25 µm |
| Vortex topological charge | 1 | 1 |
| Measured probe exit FWHM | ≈0.96 µm | ≈1.2 µm |
| Pump diameter | **0.91 µm hole width** | ≈6 µm ring diameter |

Asaf's example `asaf.m` uses a pump waist location of 17 µm and probe waist location of 20 µm. `mainProjectSimulation.m` uses a pump/probe waist location of 25 µm. The sweep starts from these code values and exposes them for modification.

### Numerical and physical validation checklist

1. Compare pump-off FWHM against the published approximately 2.77–2.8 µm for the narrow probe.
2. Check grid convergence of the <1 µm focused spot by halving radial spacing; repeat with finer axial and temporal grids.
3. Check radial-domain convergence, especially for the **100 µm plane-wave approximation**. A 30 µm domain is not sufficiently wide to faithfully represent a 100 µm Gaussian spot; use a larger domain or a true uniform entrance-field option for publication.
4. Verify that the simulated pump peak-to-peak diameter at the sample rear equals the requested x-coordinate.
5. Validate total pump energy at the entrance. The original `laser.computeE0FromEnergy` uses the Gaussian power integral even for an LG vortex, so `pump_width_sweep` applies a vortex normalization correction to the field.
6. Verify experimental pump **peak-to-peak** diameters from measured line profiles before asserting correlation.
7. Test carrier-density limitations, increased absorption, damage threshold, and the possible appearance of side lobes. FWHM alone is not a complete optical resolution metric.
8. Check Fresnel interface transmission and definitions of reported pulse energy before making absolute-parameter comparisons.

**Status:** The extension is prepared for running in MATLAB but has not been numerically executed or validated in this environment. A plotted agreement is not guaranteed.

## Known limitations

The scalar paraxial approximation and cylindrical symmetry exclude nonaxisymmetric aberrations and full-vector optics. The radial mask is a numerical absorber, not a formally matched electromagnetic PML. Carrier generation is assumed to have unity quantum efficiency and uses a phenomenological generation-time envelope. Strongly pumped silicon may require nonlinear carrier-rate equations and density-dependent optical constants. Energy conversion and idealized vortex shape should be checked against measured sample-plane profiles.

## References

1. A. Malran, D. Adar, Z. Zalevsky and M. Sinvani, *Pump-Width-Controlled Super-Resolution Focusing in Silicon with a Single-Wavelength Vortex Pump-Probe Scheme*, draft manuscript, especially Figures 3, 6 and 7.
2. N. Shabairou, Z. Zalevsky and M. Sinvani, “Laser Beam Self-Focusing in Silicon at an Absorbed Wavelength by a Vortex Beam in the Same Wavelength,” *ACS Omega* **9** (2024), 16969–16975. https://doi.org/10.1021/acsomega.3c08325
3. N. Shabairou, M. Tiferet, Z. Zalevsky and M. Sinvani, “Dynamics of laser-induced tunable focusing in silicon,” *Scientific Reports* **12** (2022), 6342. https://doi.org/10.1038/s41598-022-10112-3
4. R. A. Soref and B. R. Bennett, “Electrooptical effects in silicon,” *IEEE Journal of Quantum Electronics* **23** (1987), 123–129.
5. L. Allen et al., “Orbital angular momentum of light and the transformation of Laguerre–Gaussian laser modes,” *Physical Review A* **45** (1992), 8185–8189. https://doi.org/10.1103/PhysRevA.45.8185

## License

No explicit license was located in the repository at the time this README was drafted. If this code is intended for reuse, add a `LICENSE` file specifying the preferred license. The absence of a license does not itself authorize third-party redistribution.
