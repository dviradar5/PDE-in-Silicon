%% DOCUMENT

function [E_rz, I_rz] = propagationBPM_rz(E0_r, r, z, probe, n_complex_rz)
%propagationBPM_rz  Axisymmetric (m=0) BPM propagation in cylindrical r-z.
%
% Uses split-step:
%   (1) half-step diffraction in r with Crank–Nicolson on axisymmetric Laplacian
%   (2) full-step medium (phase + absorption) from n_complex(r,z)
%   (3) half-step diffraction
%
% Inputs:
%   E0_r         : Nr×1 complex field envelope at z=0 (vs r)
%   r            : Nr×1 radial grid [m], includes r=0, increasing
%   z            : Nz×1 propagation grid [m]
%   probe        : object/struct with probe.lambda (used for k0)
%   n_complex_rz : Nr×Nz complex refractive index map n + i*kext  (fixed time)
%
% Outputs:
%   E_rz : Nr×Nz complex field envelope
%   I_rz : Nr×Nz intensity |E|^2
%   rVec : Nr×1
%   zVec : Nz×1

    sp = systemParameters();

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

    % Uniform steps (recommended)
    dr = mean(diff(rVec));
    if any(abs(diff(rVec) - dr) > 1e-12*max(1,abs(dr)))
        warning("r grid not uniform. Using mean(dr)=%.3g m", dr);
    end

    dz = mean(diff(zVec));
    if any(abs(diff(zVec) - dz) > 1e-12*max(1,abs(dz)))
        warning("z grid not uniform. Using mean(dz)=%.3g m", dz);
    end

    k0 = 2*pi/probe.lambda;   % vacuum wavenumber
    n_ref = sp.n;             % reference real index
    k = k0 * n_ref;           % reference propagation constant

    % Build axisymmetric radial Laplacian operator L such that:
    %   (L E)_i ≈ d2E/dr2 + (1/r) dE/dr   (m=0)
    %
    % CN half-step for diffraction:
    %   (I - a L) E^{*} = (I + a L) E
    % where a = i dz/(4k)
    a = 1i * dz / (4*k);

    L = radialLaplacianMatrix(rVec, dr);      % Nr×Nr sparse
    Aimp = speye(Nr) - a*L;
    Aexp = speye(Nr) + a*L;

    % Pre-factorize for speed
    % (MATLAB will pick a sparse solver)
    % You can also use decomposition(Aimp,'lu') in newer MATLAB
    % but backslash is fine here.

    E_rz = zeros(Nr, Nz);
    E_rz(:,1) = E0_r(:);

    E = E0_r(:);

    for iz = 2:Nz
        % ---- half diffraction (CN) ----
        E = Aimp \ (Aexp * E);

        % ---- medium step (phase + absorption) ----
        n_here = n_complex_rz(:,iz);
        n_real = real(n_here);
        n_imag = imag(n_here);

        % phase from Δn relative to n_ref
        E = E .* exp( 1i * k0 * (n_real - n_ref) * dz );

        % absorption from extinction coefficient kext
        E = E .* exp( -k0 * n_imag * dz );

        % ---- half diffraction (CN) ----
        E = Aimp \ (Aexp * E);

        E_rz(:,iz) = E;
    end

    E_rz = probe.updateProfile(E_rz);
    I_rz = 0.5 * sp.n * sp.eps0 * sp.c0 * abs(E_rz).^2;

end


function L = radialLaplacianMatrix(rVec, dr)
%radialLaplacianMatrix  Axisymmetric Laplacian in r for m=0.
% Handles r=0 with symmetry: dE/dr|_{r=0} = 0

    Nr = numel(rVec);
    if abs(rVec(1)) > 1e-15
        warning("r(1) is not ~0. For best accuracy, include r=0 as first sample.");
    end

    main = zeros(Nr,1);
    up   = zeros(Nr-1,1);
    low  = zeros(Nr-1,1);

    % Interior points i=2..Nr-1
    for i = 2:Nr-1
        ri = rVec(i);
        % Standard centered FD:
        % d2/dr2: (E_{i+1}-2E_i+E_{i-1})/dr^2
        % (1/r)d/dr: (1/ri)*(E_{i+1}-E_{i-1})/(2dr)
        low(i-1) = 1/dr^2 - 1/(2*ri*dr);
        main(i)  = -2/dr^2;
        up(i)    = 1/dr^2 + 1/(2*ri*dr);
    end

    % Boundary at r=0 (i=1), Neumann: dE/dr = 0
    % Use ghost point E_0 = E_2, then:
    % d2E/dr2 at r=0 ≈ (E2 - 2E1 + E0)/dr^2 = 2(E2 - E1)/dr^2
    % and (1/r)dE/dr term -> 0 at r=0 for symmetry
    main(1) = -2/dr^2;
    up(1)   = 2/dr^2;

    % Outer boundary at r=rmax (i=Nr)
    % Simple Neumann dE/dr=0: E_{Nr+1}=E_{Nr-1}
    % => d2E/dr2 ≈ 2(E_{Nr-1}-E_{Nr})/dr^2, (1/r)dE/dr ≈ 0
    main(Nr)    = -2/dr^2;
    low(Nr-1)   = 2/dr^2;

    L = spdiags([ [low;0], main, [0;up] ], [-1,0,1], Nr, Nr);
end
