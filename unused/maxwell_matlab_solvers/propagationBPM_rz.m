function [E_rz, I_rz, info] = propagationBPM_rz(E0_r, r, z, probe, n_complex, a, b)
%PROPAGATIONBPM_RZ Axisymmetric finite-difference BPM solver.
%
% Syntax
%   [E_rz, I_rz] = propagationBPM_rz(E0_r, r, z, probe, n_complex, a, b)
%   [E_rz, I_rz, info] = propagationBPM_rz(...)
%
% Inputs
%   E0_r      : complex field envelope at z(1), size Nr x 1.
%   r         : radial grid [m], size Nr x 1. Prefer r(1)=0 and uniform dr.
%   z         : propagation grid [m], size Nz x 1.
%   probe     : struct/object containing wavelength in meters. Supported
%               names: probe.lambda, probe.wl, probe.wl0, probe.wavelength,
%               probe.wl2. Optional: probe.n0 for reference index.
%   n_complex : complex refractive index n(r,z)=n_r+i*k. Accepted sizes:
%               scalar, Nr-vector, Nz-vector, Nr x Nz, or Nz x Nr.
%   a,b       : radial sponge mask parameters: exp(-a*(r/rmax)^b).
%
% Outputs
%   E_rz      : complex field envelope, size Nr x Nz.
%   I_rz      : intensity in arbitrary units, abs(E_rz).^2.
%   info      : numerical/meta data.
%
% Equation solved
%   Write E_phys(r,z) = E_rz(r,z)*exp(i*k0*n0*z). Under the scalar,
%   forward, paraxial approximation,
%
%       dE/dz = i/(2*k0*n0) * nabla_perp^2 E
%              + i*k0/(2*n0) * (n(r,z)^2 - n0^2) * E.
%
%   The radial Laplacian is d2/dr2 + (1/r)d/dr. The z step is Crank-Nicolson.
%
% Important limits
%   Good for slowly varying envelopes, mostly forward propagation, weak
%   reflections, and paraxial angles. For high NA, sharp interfaces, or
%   back-reflections, compare to FDFD/FDTD.

    arguments
        E0_r (:,1) {mustBeNumeric}
        r (:,1) double {mustBeFinite}
        z (:,1) double {mustBeFinite}
        probe
        n_complex {mustBeNumeric}
        a double = 0
        b double = 8
    end

    Nr = numel(r);
    Nz = numel(z);
    if numel(E0_r) ~= Nr
        error('E0_r must have the same number of elements as r.');
    end
    if Nz < 2
        error('z must contain at least two points.');
    end
    if any(diff(r) <= 0) || any(diff(z) <= 0)
        error('r and z must be strictly increasing.');
    end

    lambda = getProbeWavelength(probe);
    k0 = 2*pi/lambda;
    n_rz = normalizeIndexMap_rz(n_complex, Nr, Nz);

    if (isstruct(probe) && isfield(probe,'n0')) || (isobject(probe) && isprop(probe,'n0'))
        n0 = double(probe.n0);
    else
        n0 = median(real(n_rz(:)), 'omitnan');
    end
    if ~isscalar(n0) || ~isfinite(n0) || n0 <= 0
        error('Reference index n0 must be a positive real scalar.');
    end

    L = radialLaplacianMatrix(r);
    eyeN = speye(Nr);
    mask = radialAbsorbingMask(r, a, b);

    E_rz = complex(zeros(Nr, Nz));
    E_rz(:,1) = E0_r(:);

    for iz = 1:Nz-1
        dz = z(iz+1) - z(iz);

        % Use the midpoint refractive-index perturbation in this z interval.
        n_mid = 0.5*(n_rz(:,iz) + n_rz(:,iz+1));
        V = 1i*k0/(2*n0) * (n_mid.^2 - n0^2);

        Aop = 1i/(2*k0*n0)*L + spdiags(V, 0, Nr, Nr);

        % Crank-Nicolson step: (I - dz/2 A)E_{z+dz}=(I + dz/2 A)E_z.
        lhs = eyeN - 0.5*dz*Aop;
        rhs = (eyeN + 0.5*dz*Aop) * E_rz(:,iz);
        E_next = lhs \ rhs;

        % Sponge suppresses artificial radial-boundary reflection.
        E_rz(:,iz+1) = E_next .* mask;
    end

    I_rz = abs(E_rz).^2;

    info = struct();
    info.method = 'axisymmetric finite-difference BPM, Crank-Nicolson';
    info.lambda = lambda;
    info.k0 = k0;
    info.n0 = n0;
    info.mask = mask;
    info.assumptions = {'scalar field','forward propagation','paraxial/small-angle','weak reflections'};
end
