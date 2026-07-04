function [E_rz, I_rz, info] = propagationSSFM2D_rz(E0_r, r, z, probe, n_complex, a, b)
%PROPAGATIONSSFM2D_RZ 2D angular-spectrum split-step propagation.
%
% Syntax
%   [E_rz, I_rz] = propagationSSFM2D_rz(E0_r, r, z, probe, n_complex, a, b)
%
% What this does
%   Converts the cylindrically symmetric input E(r) and n(r,z) to a 2D
%   Cartesian grid E(x,y), propagates each z step with an exact homogeneous
%   angular-spectrum operator, and applies the spatially varying refractive
%   index as split-step phase/absorption factors.
%
% Equation model
%   Scalar one-way Helmholtz split-step:
%       homogeneous diffraction: exp(i*kz*dz), kz=sqrt((k0*n0)^2-kx^2-ky^2)
%       index perturbation:      exp(i*k0*(n(x,y,z)-n0)*dz)
%
% Why use it
%   It is less paraxial than BPM for propagation through a nearly homogeneous
%   reference medium, because the homogeneous diffraction step uses exact kz.
%   It is still scalar and one-way, so it does not model polarization coupling
%   or back-reflections.

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
    if abs(r(1)) > 1e-12 * max(r(end), 1)
        error('This 2D reconstruction expects r(1)=0.');
    end

    drs = diff(r);
    dr = mean(drs);
    if max(abs(drs - dr)) > 1e-9 * max(abs(dr), eps)
        error('This solver expects uniform r spacing.');
    end

    lambda = getProbeWavelength(probe);
    k0 = 2*pi/lambda;
    n_rz = normalizeIndexMap_rz(n_complex, Nr, Nz);

    if (isstruct(probe) && isfield(probe,'n0')) || (isobject(probe) && isprop(probe,'n0'))
        n0 = double(probe.n0);
    else
        n0 = median(real(n_rz(:)), 'omitnan');
    end

    % Symmetric Cartesian grid: x = [-rmax ... 0 ... rmax].
    x = [-flipud(r(2:end)); r];
    Nx = numel(x);
    [X, Y] = meshgrid(x, x);
    RHO = hypot(X, Y);
    mid = Nr;  % index of x=0 and y=0 in the symmetric grid

    E_xy = interp1(r, E0_r(:), RHO, 'linear', 0);
    n_edge = n_rz(end,:);

    % FFT-space coordinates in MATLAB FFT ordering.
    dk = 2*pi/(Nx*dr);
    if mod(Nx,2) == 0
        kvals = [0:(Nx/2-1), -Nx/2:-1] * dk;
    else
        kvals = [0:((Nx-1)/2), -((Nx-1)/2):-1] * dk;
    end
    [KX, KY] = meshgrid(kvals, kvals);
    KZ = sqrt(complex((k0*n0)^2 - KX.^2 - KY.^2, 0));

    % 2D radial sponge. Stronger outside the original radial disk.
    if nargin < 6 || isempty(a), a = 0; end
    if nargin < 7 || isempty(b), b = 8; end
    if a == 0
        mask2 = ones(size(RHO));
    else
        mask2 = exp(-a * (RHO / max(r)).^b);
    end

    E_rz = complex(zeros(Nr, Nz));
    E_rz(:,1) = E0_r(:);

    for iz = 1:Nz-1
        dz = z(iz+1) - z(iz);

        n1 = interp1(r, n_rz(:,iz),   RHO, 'linear', n_edge(iz));
        n2 = interp1(r, n_rz(:,iz+1), RHO, 'linear', n_edge(iz+1));

        % Symmetric split step: half index, homogeneous diffraction, half index.
        E_xy = exp(1i*k0*(n1 - n0)*dz/2) .* E_xy;
        E_xy = ifft2(fft2(E_xy) .* exp(1i*KZ*dz));
        E_xy = exp(1i*k0*(n2 - n0)*dz/2) .* E_xy;

        E_xy = E_xy .* mask2;

        % Extract the x>=0 radial line at y=0.
        E_rz(:,iz+1) = E_xy(mid, mid:end).';
    end

    I_rz = abs(E_rz).^2;

    info = struct();
    info.method = '2D angular-spectrum split-step scalar propagation';
    info.lambda = lambda;
    info.k0 = k0;
    info.n0 = n0;
    info.x = x;
    info.assumptions = {'scalar field','one-way propagation','radially symmetric input/index','weak back-reflection'};
end
