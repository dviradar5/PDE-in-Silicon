%% FINISHED

function [E_rz, I_rz] = maxwell_fdfd_corrected(E0_r, r, z, probe, n_complex, opts)
    % Scalar FDFD propagation calculator
    % ---------------------------------------------------------------------
    % Calculates the probe's field inside the sample by solving the scalar
    % frequency-domain Helmholtz equation in cylindrical coordinates:
    %
    %       (1/r*d/dr(r*dE/dr) + d^2E/dz^2 + k0^2*n(r,z)^2*E = 0)
    %
    % This is meant to be compared with propagationBPM_rz.  The inputs and
    % outputs are intentionally the same as propagationBPM_rz:
    %
    % =====================================================================
    % INPUTS:
    %        E0_r - complex electric field vector at z=0, Nr [V/m]
    %        r - radial coordinate vector, Nr [m]
    %        z - z coordinate vector, Nz [m]
    %        probe - probe laser beam, Laser-type object
    %        n_complex - complex refractive index matrix, Nr x Nz
    %                    convention: n_complex = n + 1i*k, k >= 0 for loss
    %        opts - optional structure:
    %               opts.npml_r       - extra radial PML grid points
    %               opts.npml_z       - extra z PML grid points after z(end)
    %               opts.pml_kmax     - maximum added extinction coefficient
    %               opts.pml_order    - PML polynomial order
    %               opts.verbose      - true/false
    %
    % OUTPUTS:
    %        E_rz - complex field inside the original sample grid, Nr x Nz [V/m]
    %        I_rz - intensity inside the original sample grid, Nr x Nz [W/m^2]
    % *********************************************************************

    sp = systemParameters();

    if nargin < 6 || isempty(opts)
        opts = struct();
    end

    rVec = r(:);
    zVec = z(:);
    E0_r = E0_r(:);

    [Nr, Nz] = size(n_complex);

    if numel(rVec) ~= Nr
        error('Number of r points must match size(n_complex,1).');
    end

    if numel(zVec) ~= Nz
        error('Number of z points must match size(n_complex,2).');
    end

    if numel(E0_r) ~= Nr
        error('E0_r must have the same number of points as r.');
    end

    % Steps:
    dr = rVec(2)-rVec(1);
    dz = zVec(2)-zVec(1);

    if max(abs(diff(rVec)-dr)) > 1e-12*abs(dr)
        error('r must be uniformly spaced.');
    end

    if max(abs(diff(zVec)-dz)) > 1e-12*abs(dz)
        error('z must be uniformly spaced.');
    end

    if abs(rVec(1)) > 1e-15
        error('For cylindrical FDFD, r(1) must be 0.');
    end

    % Defaults:
    if ~isfield(opts,'npml_r')
        opts.npml_r = max(10, round(0.15*Nr));
    end

    if ~isfield(opts,'npml_z')
        opts.npml_z = max(10, round(0.15*Nz));
    end

    if ~isfield(opts,'pml_kmax')
        opts.pml_kmax = 1.0;
    end

    if ~isfield(opts,'pml_order')
        opts.pml_order = 3;
    end

    if ~isfield(opts,'verbose')
        opts.verbose = true;
    end

    npml_r = opts.npml_r;
    npml_z = opts.npml_z;

    % ------------------------------------------------------------------
    % Extend the computational domain.
    %
    % Important:
    % The PML is added OUTSIDE the original r,z grid, not inside it.
    % Therefore the returned E_rz and I_rz are not directly damped by the PML.
    % ------------------------------------------------------------------
    NrT = Nr + npml_r;
    NzT = Nz + npml_z;

    rExt = (0:NrT-1).' * dr;
    zExt = zVec(1) + (0:NzT-1).' * dz;

    nExt = zeros(NrT, NzT);

    % Original physical region:
    nExt(1:Nr,1:Nz) = n_complex;

    % Radial extension: continue the last physical radial value:
    if npml_r > 0
        for iz = 1:Nz
            nExt(Nr+1:NrT,iz) = n_complex(end,iz);
        end
    end

    % z extension: continue the last physical z value:
    if npml_z > 0
        for ir = 1:NrT
            nExt(ir,Nz+1:NzT) = nExt(ir,Nz);
        end
    end

    % Add lossy PML only in the added extension region:
    nExt = add_pml_loss(nExt, Nr, Nz, npml_r, npml_z, ...
                        opts.pml_kmax, opts.pml_order);

    % Boundary field at z = z(1):
    E0_ext = zeros(NrT,1);
    E0_ext(1:Nr) = E0_r;

    k0 = 2*pi/probe.lambda;

    if opts.verbose
        fprintf('FDFD grid: Nr x Nz = %d x %d, total unknowns = %d\n', ...
                NrT, NzT, NrT*NzT);
    end

    % Flattening convention:
    % p = ir + (iz-1)*NrT
    ind = @(ir,iz) ir + (iz-1)*NrT;

    Ntot = NrT*NzT;

    % Preallocate approximately 5 entries per row:
    rows = zeros(5*Ntot,1);
    cols = zeros(5*Ntot,1);
    vals = zeros(5*Ntot,1);
    count = 0;

    b = zeros(Ntot,1);

    inv2rdr = zeros(NrT,1);
    inv2rdr(2:end) = 1 ./ (2*rExt(2:end)*dr);

    ar = 1/dr^2;
    az = 1/dz^2;

    for iz = 1:NzT
        for ir = 1:NrT

            p = ind(ir,iz);

            % ----------------------------------------------------------
            % Input boundary:
            % force E(r,z_start) = E0(r).
            % This is the FDFD equivalent of launching the BPM initial field.
            % ----------------------------------------------------------
            if iz == 1
                count = count + 1;
                rows(count) = p;
                cols(count) = p;
                vals(count) = 1;
                b(p) = E0_ext(ir);
                continue;
            end

            % ----------------------------------------------------------
            % Outer boundaries:
            % Dirichlet at the OUTER edge of the extended PML.
            % The field should already be attenuated before reaching it.
            % ----------------------------------------------------------
            if ir == NrT || iz == NzT
                count = count + 1;
                rows(count) = p;
                cols(count) = p;
                vals(count) = 1;
                b(p) = 0;
                continue;
            end

            n2 = nExt(ir,iz)^2;

            % ----------------------------------------------------------
            % Cylindrical axis r=0:
            %
            % For an axisymmetric field:
            %   lap_r(E)|_{r=0} = d2E/dr2 + (1/r)dE/dr = 4*(E2-E1)/dr^2
            %
            % This is a common source of factor-of-2 errors.
            % ----------------------------------------------------------
            if ir == 1
                rr = 4/dr^2;

                % Center
                count = count + 1;
                rows(count) = p;
                cols(count) = p;
                vals(count) = -rr - 2*az + k0^2*n2;

                % r neighbor
                count = count + 1;
                rows(count) = p;
                cols(count) = ind(ir+1,iz);
                vals(count) = rr;

                % z backward
                count = count + 1;
                rows(count) = p;
                cols(count) = ind(ir,iz-1);
                vals(count) = az;

                % z forward
                count = count + 1;
                rows(count) = p;
                cols(count) = ind(ir,iz+1);
                vals(count) = az;

                continue;
            end

            % ----------------------------------------------------------
            % Interior cylindrical Helmholtz stencil:
            %
            %   d2E/dr2 + (1/r)dE/dr + d2E/dz2 + k0^2*n^2*E = 0
            % ----------------------------------------------------------
            a_plus  = ar + inv2rdr(ir);
            a_minus = ar - inv2rdr(ir);
            a_center = -2*ar - 2*az + k0^2*n2;

            % Center
            count = count + 1;
            rows(count) = p;
            cols(count) = p;
            vals(count) = a_center;

            % r + dr
            count = count + 1;
            rows(count) = p;
            cols(count) = ind(ir+1,iz);
            vals(count) = a_plus;

            % r - dr
            count = count + 1;
            rows(count) = p;
            cols(count) = ind(ir-1,iz);
            vals(count) = a_minus;

            % z + dz
            count = count + 1;
            rows(count) = p;
            cols(count) = ind(ir,iz+1);
            vals(count) = az;

            % z - dz
            count = count + 1;
            rows(count) = p;
            cols(count) = ind(ir,iz-1);
            vals(count) = az;
        end
    end

    rows = rows(1:count);
    cols = cols(1:count);
    vals = vals(1:count);

    A = sparse(rows, cols, vals, Ntot, Ntot);

    if opts.verbose
        fprintf('Solving scalar FDFD Helmholtz system...\n');
    end

    Evec = A \ b;

    E_ext = reshape(Evec, NrT, NzT);

    % Return only the original physical region:
    E_rz = E_ext(1:Nr,1:Nz);

    % Match propagationBPM_rz intensity convention:
    I_rz = 0.5 * sp.n * sp.eps0 * sp.c0 * abs(E_rz).^2;

end


function nExt = add_pml_loss(nExt, Nr, Nz, npml_r, npml_z, pml_kmax, pml_order)
    % Adds artificial extinction coefficient outside the physical region.
    % n_complex convention: n = n_real + 1i*k, k > 0 gives absorption.

    [NrT, NzT] = size(nExt);

    if npml_r > 0
        for ir = Nr+1:NrT
            s = (ir - Nr) / npml_r;
            add_k = pml_kmax * s^pml_order;
            nExt(ir,:) = nExt(ir,:) + 1i*add_k;
        end
    end

    if npml_z > 0
        for iz = Nz+1:NzT
            s = (iz - Nz) / npml_z;
            add_k = pml_kmax * s^pml_order;
            nExt(:,iz) = nExt(:,iz) + 1i*add_k;
        end
    end
end
