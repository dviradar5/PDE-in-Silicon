%% DOCUMENT

function [E_xz, I_xz, xVec, zVec] = propagationBPM(E0_r, r, z, probe, n_complex_rz, useInterp)
%BPM_probe_silicon_xz  Split-step Fourier BPM for probe propagation in silicon (x-z).
%
% Inputs:
%   E0_r          : Nr×1 complex field envelope at z=0 in cylindrical r
%   r             : Nr×1 radial grid (>=0, increasing)
%   z             : Nz×1 propagation grid
%   lambda        : wavelength [m]
%   n_complex_rz  : Nr×Nz complex refractive index map (n + i*k) at a fixed time
%   sp.n         : reference real refractive index for BPM (e.g. sp.n)
%   useInterp     : true/false. If false tries exact mapping r<->|x|; else interp1
%
% Outputs:
%   E_xz : Nx×Nz complex field envelope in Cartesian x-z
%   I_xz : Nx×Nz intensity |E|^2
%   xVec : Nx×1 x grid (symmetric from r)
%   zVec : Nz×1 z grid
    
    sp = systemParameters();
    
    if nargin < 5 || isempty(sp.n), sp.n = real(n_complex_rz(1,1)); end
    if nargin < 6 || isempty(useInterp), useInterp = true; end

    rVec = r(:);
    zVec = z(:);

    [Nr, Nz] = size(n_complex_rz);
    if numel(rVec) ~= Nr
        error("length(r) (%d) must match size(n_complex_rz,1) (%d)", numel(rVec), Nr);
    end
    if numel(zVec) ~= Nz
        error("length(z) (%d) must match size(n_complex_rz,2) (%d)", numel(zVec), Nz);
    end
    if numel(E0_r) ~= Nr
        error("E0_r must be Nr×1 (Nr=%d)", Nr);
    end

    % ----- build symmetric x grid from r -----
    xVec = [-flipud(rVec(2:end)); rVec];    % length Nx = 2*Nr-1, x=0 at center
    absx = abs(xVec);
    Nx = numel(xVec);

    % ----- map E0(r) -> E0(x) using symmetry -----
    if useInterp
        E0_x = interp1(rVec, E0_r(:), absx, 'linear', 0);
    else
        tol = 1e-12 * max(1, max(rVec));
        [tf, idx] = ismembertol(absx, rVec, tol);
        if ~all(tf)
            E0_x = interp1(rVec, E0_r(:), absx, 'linear', 0);
        else
            E0_x = E0_r(idx);
        end
    end

    % ----- map n(r,z) -> n(x,z) (complex) -----
    n_xz = zeros(Nx, Nz);
    if useInterp
        for iz = 1:Nz
            n_xz(:,iz) = interp1(rVec, n_complex_rz(:,iz), absx, 'linear', 'extrap');
        end
    else
        tol = 1e-12 * max(1, max(rVec));
        [tf, idx] = ismembertol(absx, rVec, tol);
        if ~all(tf)
            for iz = 1:Nz
                n_xz(:,iz) = interp1(rVec, n_complex_rz(:,iz), absx, 'linear', 'extrap');
            end
        else
            n_xz = n_complex_rz(idx, :);
        end
    end

    % ----- BPM setup -----
    % assume (approximately) uniform dz
    dz = mean(diff(zVec));
    if any(abs(diff(zVec) - dz) > 1e-12*max(1,abs(dz)))
        warning("z grid not uniform. Using mean(dz)=%.3g m", dz);
    end

    dx = mean(diff(xVec));
    k0 = 2*pi/probe.lambda;              % vacuum wavenumber
    k  = k0 * sp.n;               % reference propagation constant in medium

    % spatial frequency grid (FFT)
    fx = (-floor(Nx/2):ceil(Nx/2)-1) / (Nx*dx);  % cycles/m
    kx = 2*pi*fx;                                % rad/m
    kx = fftshift(kx).';                         % column, aligned with fft

    % diffraction operator (half-step)
    Hhalf = exp(-1i * (kx.^2) * (dz/(2*k)));      % exp[-i kx^2 dz/(2k)]

    % output arrays
    E_xz = zeros(Nx, Nz);
    E_xz(:,1) = E0_x;

    % ----- split-step propagation -----
    A = E0_x;

    for iz = 2:Nz
        % half diffraction
        A = ifft( fft(A) .* Hhalf );

        % medium (phase + absorption) using n(x,z) at current plane
        n_here = n_xz(:,iz);

        % n = n' + i*kext
        n_real = real(n_here);
        n_imag = imag(n_here);

        % phase from index perturbation relative to sp.n
        A = A .* exp( 1i * k0 * (n_real - sp.n) * dz );

        % absorption from extinction coefficient kext: exp(-k0*kext*dz)
        A = A .* exp( -k0 * n_imag * dz );

        % half diffraction
        A = ifft( fft(A) .* Hhalf );

        E_xz(:,iz) = A;
    end
    
    probe.updateProfile(E_xz)
    I_xz = probe.intensityProfileBLDumped(z);
end
