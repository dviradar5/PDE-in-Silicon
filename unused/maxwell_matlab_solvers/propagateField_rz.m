function [E_rz, I_rz, info] = propagateField_rz(method, E0_r, r, z, probe, n_complex, a, b)
%PROPAGATEFIELD_RZ Dispatcher for comparing propagation solvers.
%
% Usage
%   [E,I] = propagateField_rz('bpm',  E0_r,r,z,probe,n_complex,a,b);
%   [E,I] = propagateField_rz('ssfm', E0_r,r,z,probe,n_complex,a,b);
%   [E,I] = propagateField_rz('fdfd', E0_r,r,z,probe,n_complex,a,b);
%
% Methods
%   'bpm'  : axisymmetric finite-difference paraxial BPM.
%   'ssfm' : 2D angular-spectrum split-step, scalar one-way.
%   'fdfd' : scalar 2D Helmholtz FDFD comparison solver.

    if nargin < 7, a = 0; end
    if nargin < 8, b = 8; end

    switch lower(string(method))
        case {"bpm", "fd-bpm", "paraxial"}
            [E_rz, I_rz, info] = propagationBPM_rz(E0_r, r, z, probe, n_complex, a, b);
        case {"ssfm", "asm", "angular", "angular-spectrum"}
            [E_rz, I_rz, info] = propagationSSFM2D_rz(E0_r, r, z, probe, n_complex, a, b);
        case {"fdfd", "helmholtz"}
            [E_rz, I_rz, info] = propagationHelmholtzFDFD_rz(E0_r, r, z, probe, n_complex, a, b);
        otherwise
            error('Unknown method "%s". Use bpm, ssfm, or fdfd.', method);
    end
end
