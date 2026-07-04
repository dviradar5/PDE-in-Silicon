# MATLAB propagation solvers for `E_rz, I_rz`

This folder contains independent MATLAB functions with the same basic inputs as:

```matlab
[E_rz, I_rz] = propagationBPM_rz(E0_r, r, z, probe, n_complex, a, b)
```

## Files

- `propagationBPM_rz.m` — axisymmetric finite-difference BPM, Crank-Nicolson in z.
- `propagationSSFM2D_rz.m` — scalar 2D angular-spectrum split-step method.
- `propagationHelmholtzFDFD_rz.m` — scalar 2D Helmholtz/FDFD comparison solver.
- `propagateField_rz.m` — dispatcher: `'bpm'`, `'ssfm'`, `'fdfd'`.
- `demo_compare_solvers.m` — minimal test script.
- Helpers: `getProbeWavelength.m`, `normalizeIndexMap_rz.m`, `radialLaplacianMatrix.m`, `radialAbsorbingMask.m`.

## Important interpretation

`I_rz = abs(E_rz).^2` is returned in arbitrary units, matching the common numerical-optics convention. If your electric-field phasor is physically calibrated in V/m, use:

```matlab
I_phys = 0.5 * real(n_complex) * eps0 * c0 .* abs(E_rz).^2;
```

## Which solver should you use?

- Use BPM for fast pump-probe parameter scans and mostly forward/paraxial propagation.
- Use SSFM/angular spectrum to check whether BPM is failing because the effective NA/angles are too large.
- Use FDFD only on small grids as a sanity check for non-paraxial/back-reflection effects.
- For real vector Maxwell in your project, use your existing `maxwell_run_Loop` / MaxwellFDFD-based route.
