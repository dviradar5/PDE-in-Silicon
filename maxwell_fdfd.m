function [E_rz, I_rz] = maxwell_fdfd(r, z, lambda, n_complex_rz, E0_r, opts)
% maxwell_fdfd  Scalar frequency-domain Helmholtz FDFD in (r,z).
%
% Solves:
%   ( 1/r d/dr (r dE/dr) + d^2E/dz^2 + k0^2 n(r,z)^2 ) E = 0
% with a "soft source" boundary at z=0: E(r, z=0) = E0_r,
% and PML-like absorption near edges (implemented as extra imag(n)).
%
% Inputs:
%   r, z           grids (uniform)
%   lambda         wavelength [m]
%   n_complex_rz   complex refractive index n + i*k  (IMPORTANT convention!)
%   E0_r           boundary field at z=0 (Nr x 1)
%   opts           struct with optional fields:
%       .npml_r    (default 20)   number of radial PML points
%       .npml_z    (default 20)   number of z PML points (applied near z=zmax)
%       .pml_kmax  (default 2.0)  max added extinction coefficient in PML (dimensionless k)
%       .bc_outer  (default "dirichlet") outer r boundary ("dirichlet" only here)
%
% Outputs:
%   E_rz, I_rz

    arguments
        r (:,1) double {mustBeNonnegative}
        z (:,1) double
        lambda (1,1) double {mustBePositive}
        n_complex_rz (:,:) double
        E0_r (:,1) double
        opts.npml_r (1,1) double {mustBeInteger, mustBeNonnegative} = 20
        opts.npml_z (1,1) double {mustBeInteger, mustBeNonnegative} = 20
        opts.pml_kmax (1,1) double {mustBeNonnegative} = 2.0
        opts.bc_outer (1,1) string = "dirichlet"
    end

    Nr = numel(r);
    Nz = numel(z);

    if size(n_complex_rz,1) ~= Nr || size(n_complex_rz,2) ~= Nz
        error("n_complex_rz must be Nr x Nz to match r and z.");
    end
    if numel(E0_r) ~= Nr
        error("E0_r must be Nr x 1 to match r.");
    end

    dr = r(2)-r(1);
    dz = z(2)-z(1);
    if any(abs(diff(r)-dr) > 1e-15), error("r must be uniformly spaced."); end
    if any(abs(diff(z)-dz) > 1e-15), error("z must be uniformly spaced."); end

    k0 = 2*pi/lambda;

    % --- Add simple "PML" by increasing extinction k near boundaries ---
    % Convention required: n_complex = n + i*k, with k >= 0 for loss.
    n_eff = add_soft_pml_to_index(r, z, n_complex_rz, opts.npml_r, opts.npml_z, opts.pml_kmax);

    % Flatten (r,z) -> linear index
    % idx(ir,iz) = ir + (iz-1)*Nr
    N = Nr*Nz;

    rows = [];
    cols = [];
    vals = [];
    b = zeros(N,1);

    % Helper inline for indexing
    idx = @(ir,iz) ir + (iz-1)*Nr;

    % Precompute 1/r terms for cylindrical Laplacian (avoid division by zero at r=0)
    inv2rdr = zeros(Nr,1);
    inv2rdr(2:end) = 1./(2*r(2:end)*dr);
    % at r=0 use symmetry: dE/dr = 0 => handled by special stencil

    for iz = 1:Nz
        for ir = 1:Nr
            p = idx(ir,iz);

            % --- z=0 boundary: enforce E(r,1) = E0_r ---
            if iz == 1
                rows(end+1) = p; cols(end+1) = p; vals(end+1) = 1;
                b(p) = E0_r(ir);
                continue;
            end

            % --- Outer boundaries: set Dirichlet E=0 (PML should make it non-reflecting-ish) ---
            if (ir == Nr) || (iz == Nz)
                rows(end+1) = p; cols(end+1) = p; vals(end+1) = 1;
                b(p) = 0;
                continue;
            end

            % --- r=0 symmetry: dE/dr = 0 ---
            if ir == 1
                % Use "mirror" point: E(-dr) = E(+dr)
                % Cyl Laplacian at r=0 becomes ~ 2*(E2 - E1)/dr^2
                % We'll implement:
                % (d2/dr2 + (1/r)d/dr)E|r=0 = 2*(E(2)-E(1))/dr^2
                rr = 2/dr^2;
                % z second derivative central
                zz = 1/dz^2;

                n2 = (n_eff(ir,iz))^2;

                % Center
                rows(end+1)=p; cols(end+1)=p;
                vals(end+1)= (-rr) + (-2*zz) + (k0^2*n2);

                % r neighbor (ir=2)
                rows(end+1)=p; cols(end+1)=idx(ir+1,iz);
                vals(end+1)= (rr);

                % z neighbors
                rows(end+1)=p; cols(end+1)=idx(ir,iz-1);
                vals(end+1)= (zz);
                rows(end+1)=p; cols(end+1)=idx(ir,iz+1);
                vals(end+1)= (zz);

                continue;
            end

            % --- Interior stencil for cylindrical Laplacian + z-Laplacian ---
            % Radial operator:
            % (1/r)d/dr(r dE/dr) ≈ (E_{i+1}-2E_i+E_{i-1})/dr^2 + (E_{i+1}-E_{i-1})/(2 r_i dr)
            a_plus  = 1/dr^2 + inv2rdr(ir);
            a_0r    = -2/dr^2;
            a_minus = 1/dr^2 - inv2rdr(ir);

            % z operator
            a_z = 1/dz^2;

            n2 = (n_eff(ir,iz))^2;

            % Center
            rows(end+1)=p; cols(end+1)=p;
            vals(end+1)= a_0r + (-2*a_z) + (k0^2*n2);

            % r neighbors
            rows(end+1)=p; cols(end+1)=idx(ir+1,iz);
            vals(end+1)= a_plus;

            rows(end+1)=p; cols(end+1)=idx(ir-1,iz);
            vals(end+1)= a_minus;

            % z neighbors
            rows(end+1)=p; cols(end+1)=idx(ir,iz+1);
            vals(end+1)= a_z;

            rows(end+1)=p; cols(end+1)=idx(ir,iz-1);
            vals(end+1)= a_z;
        end
    end

    A = sparse(rows, cols, vals, N, N);

    % Solve
    Evec = A \ b;
    E_rz = reshape(Evec, Nr, Nz);
    I_rz = abs(E_rz).^2;
end


function n_eff = add_soft_pml_to_index(r, z, n_complex_rz, npml_r, npml_z, kmax)
% Add an extra extinction coefficient near r=Rmax and z=Zmax.
% This is NOT a true coordinate-stretched PML, but a practical absorber.

    n_eff = n_complex_rz;

    Nr = numel(r);
    Nz = numel(z);

    if npml_r > 0
        ir0 = max(1, Nr - npml_r + 1);
        s = linspace(0,1, Nr-ir0+1).';   % 0..1
        prof = (s.^2);                   % quadratic ramp
        addk_r = kmax * prof;            % dimensionless extinction
        for ir = ir0:Nr
            n_eff(ir,:) = n_eff(ir,:) + 1i*addk_r(ir-ir0+1);
        end
    end

    if npml_z > 0
        iz0 = max(1, Nz - npml_z + 1);
        s = linspace(0,1, Nz-iz0+1);     % 0..1
        prof = (s.^2);
        addk_z = kmax * prof;
        for iz = iz0:Nz
            n_eff(:,iz) = n_eff(:,iz) + 1i*addk_z(iz-iz0+1);
        end
    end
end
