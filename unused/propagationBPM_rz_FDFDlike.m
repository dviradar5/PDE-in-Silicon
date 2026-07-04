function [E_rz, I_rz] = propagationBPM_rz_FDFDlike(E0_r, r, z, probe, n_complex, a, b)
    % propagationBPM_rz_FDFDlike
    % ---------------------------------------------------------------------
    % One-way paraxial BPM in cylindrical coordinates, written to be closer
    % to the scalar Helmholtz equation used in maxwell_fdfd_corrected:
    %
    %   dE/dz = i/(2*k0*n0) * Laplacian_r(E)
    %          + i*k0*(n^2-n0^2)/(2*n0) * E
    %
    % n_complex convention:
    %   n_complex = n_real + 1i*k, with k >= 0 for loss.
    %
    % INPUTS:
    %   E0_r      - input field at z(1), Nr x 1
    %   r         - radial vector [m], starts at 0
    %   z         - axial vector [m]
    %   probe     - laser object with probe.lambda
    %   n_complex - Nr x Nz complex refractive index
    %   a,b       - radial damping mask parameters
    %
    % OUTPUTS:
    %   E_rz      - propagated field
    %   I_rz      - intensity [W/m^2], using same convention as FDFD:
    %               I = 0.5*eps0*c0*n0*|E|^2

    sp = systemParameters();

    rVec = r(:);
    zVec = z(:);
    E = E0_r(:);

    [Nr, Nz] = size(n_complex);

    if numel(rVec) ~= Nr
        error("length(r) must match size(n_complex,1).");
    end

    if numel(zVec) ~= Nz
        error("length(z) must match size(n_complex,2).");
    end

    if numel(E) ~= Nr
        error("length(E0_r) must match length(r).");
    end

    if abs(rVec(1)) > 1e-15
        error("For cylindrical BPM, r(1) must be 0.");
    end

    if any(imag(n_complex(:)) < -1e-12)
        warning("imag(n_complex)<0 found. This BPM assumes n = n_real + 1i*k, k >= 0 for loss.");
    end

    dr = rVec(2) - rVec(1);
    dz = zVec(2) - zVec(1);

    k0 = 2*pi/probe.lambda;
    n0 = sp.n;   % reference/background silicon index

    % Cylindrical radial Laplacian using the same r=0 factor as FDFD:
    L = lapCylR(rVec, dr);

    % True half-step diffraction CN coefficient:
    sHalf = 1i * dz / (8*n0*k0);

    Aimp = speye(Nr) - sHalf * L;
    Aexp = speye(Nr) + sHalf * L;

    % Radial damping mask, only near large r:
    if nargin < 6 || isempty(a)
        a = 0;
    end
    if nargin < 7 || isempty(b)
        b = 8;
    end
    mask = exp(-a * (rVec / rVec(end)).^b);

    E_rz = complex(zeros(Nr, Nz));
    I_rz = zeros(Nr, Nz);

    E_rz(:,1) = E;

    % Match maxwell_fdfd_corrected intensity convention:
    I_rz(:,1) = 0.5 * sp.eps0 * sp.c0 * n0 .* abs(E).^2;

    for iz = 2:Nz

        % Use midpoint refractive index in z.
        n_mid = 0.5 * (n_complex(:,iz-1) + n_complex(:,iz));

        % Half-step diffraction:
        E = Aimp \ (Aexp * E);

        % Full medium step from scalar Helmholtz envelope equation.
        % This includes both phase and absorption.
        deltaBeta = k0 * (n_mid.^2 - n0^2) / (2*n0);

        E = E .* exp(1i * deltaBeta * dz);

        % Half-step diffraction:
        E = Aimp \ (Aexp * E);

        % Boundary damping:
        E = E .* mask;

        E_rz(:,iz) = E;

        % Use same intensity convention as FDFD for direct comparison.
        I_rz(:,iz) = 0.5 * sp.eps0 * sp.c0 * n0 .* abs(E).^2;
    end
end


function L = lapCylR(rVec, dr)
    % Cylindrical radial Laplacian:
    %
    %   L(E) = d2E/dr2 + (1/r)dE/dr
    %
    % At r=0:
    %   L(E)|_0 = 4*(E2-E1)/dr^2
    %
    % This matches the axis treatment in maxwell_fdfd_corrected.

    Nr = numel(rVec);

    rows = [];
    cols = [];
    vals = [];

    add = @(i,j,v) assignin("caller", "dummy", 0);

    % Helper for appending sparse entries:
    function push(i,j,v)
        rows(end+1,1) = i;
        cols(end+1,1) = j;
        vals(end+1,1) = v;
    end

    % r = 0 axis:
    push(1,1, -4/dr^2);
    push(1,2,  4/dr^2);

    % Interior:
    for ir = 2:Nr-1
        ri = rVec(ir);

        a_plus  = 1/dr^2 + 1/(2*ri*dr);
        a_minus = 1/dr^2 - 1/(2*ri*dr);
        a0      = -2/dr^2;

        push(ir, ir-1, a_minus);
        push(ir, ir,   a0);
        push(ir, ir+1, a_plus);
    end

    % Outer boundary: simple Neumann dE/dr = 0.
    % The external radial mask should already suppress the field near here.
    push(Nr, Nr-1,  2/dr^2);
    push(Nr, Nr,   -2/dr^2);

    L = sparse(rows, cols, vals, Nr, Nr);
end