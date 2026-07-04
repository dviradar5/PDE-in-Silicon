function [E_rz, I_rz, info] = propagationHelmholtzFDFD_rz(E0_r, r, z, probe, n_complex, a, b)
%PROPAGATIONHELMHOLTZFDFD_RZ Scalar 2D Helmholtz/FDFD comparison solver.
%
% Syntax
%   [E_rz, I_rz] = propagationHelmholtzFDFD_rz(E0_r, r, z, probe, n_complex, a, b)
%
% Model solved
%   For a single scalar field component, e.g. TE Ey in an x-z invariant 2D
%   model with mu=mu0 and no vector coupling:
%
%       d2E/dx2 + d2E/dz2 + k0^2*n(x,z)^2*E = 0.
%
%   This is a frequency-domain finite-difference Helmholtz equation. It is
%   closer to Maxwell than BPM because it is not a forward-paraxial marching
%   equation, but it is NOT a full vector Maxwell solver.
%
% Boundary conditions
%   z=z(1):    Dirichlet injection E(x,z1)=E0(|x|)
%   x=edges:   Dirichlet E=0, intended to be far from the beam; optional
%              sponge loss is added near edges when a>0.
%   z=z(end):  first-order outgoing condition dE/dz = i*k0*n*E.
%
% Practical warning
%   FDFD builds and solves one sparse linear system of size Nx*Nz. Use it on
%   coarser/smaller grids for validation. For production vector Maxwell, use
%   MaxwellFDFD/your maxwell_run_Loop path.

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
        error('This FDFD x-z reconstruction expects r(1)=0.');
    end

    dr = mean(diff(r));
    if max(abs(diff(r) - dr)) > 1e-9 * max(abs(dr), eps)
        error('This FDFD implementation expects uniform r spacing.');
    end
    dzs = diff(z);
    dz = mean(dzs);
    if max(abs(dzs - dz)) > 1e-9 * max(abs(dz), eps)
        error('This FDFD implementation expects uniform z spacing.');
    end

    lambda = getProbeWavelength(probe);
    k0 = 2*pi/lambda;
    n_rz = normalizeIndexMap_rz(n_complex, Nr, Nz);

    x = [-flipud(r(2:end)); r];
    Nx = numel(x);
    E0_x = interp1(r, E0_r(:), abs(x), 'linear', 0);

    % Build n(x,z) from n(r,z).
    n_xz = complex(zeros(Nx, Nz));
    for j = 1:Nz
        n_xz(:,j) = interp1(r, n_rz(:,j), abs(x), 'linear', n_rz(end,j));
    end

    % Optional lossy sponge near transverse edges and output boundary.
    if nargin >= 6 && ~isempty(a) && a > 0
        sx = abs(x) / max(abs(x));
        loss_x = a * max(0, (sx - 0.75) / 0.25).^b;
        sz = ((1:Nz) - 1) / max(Nz-1, 1);
        loss_z = a * max(0, (sz - 0.85) / 0.15).^b;
        n_xz = n_xz + 1i*(loss_x(:) + loss_z(:).');
    end

    N = Nx * Nz;
    maxEntries = 5*N;
    rows = zeros(maxEntries,1);
    cols = zeros(maxEntries,1);
    vals = complex(zeros(maxEntries,1));
    rhs  = complex(zeros(N,1));
    nnzCount = 0;

    idx = @(i,j) i + (j-1)*Nx;

    function add(row, col, val)
        nnzCount = nnzCount + 1; %#ok<NASGU>
        rows(nnzCount) = row;
        cols(nnzCount) = col;
        vals(nnzCount) = val;
    end

    for j = 1:Nz
        for i = 1:Nx
            p = idx(i,j);

            if j == 1
                % Input plane: impose the requested initial field.
                add(p, p, 1);
                rhs(p) = E0_x(i);

            elseif i == 1 || i == Nx
                % Transverse computational edge; keep it far + sponge.
                add(p, p, 1);
                rhs(p) = 0;

            elseif j == Nz
                % First-order outgoing boundary: (E_j - E_{j-1})/dz = i*k*E_j.
                nloc = n_xz(i,j);
                add(p, p, 1 - 1i*k0*nloc*dz);
                add(p, idx(i,j-1), -1);
                rhs(p) = 0;

            else
                % Interior five-point Helmholtz stencil.
                nloc = n_xz(i,j);
                add(p, idx(i-1,j), 1/dr^2);
                add(p, idx(i+1,j), 1/dr^2);
                add(p, idx(i,j-1), 1/dz^2);
                add(p, idx(i,j+1), 1/dz^2);
                add(p, p, -2/dr^2 - 2/dz^2 + (k0*nloc)^2);
                rhs(p) = 0;
            end
        end
    end

    rows = rows(1:nnzCount);
    cols = cols(1:nnzCount);
    vals = vals(1:nnzCount);
    A = sparse(rows, cols, vals, N, N);

    E_vec = A \ rhs;
    E_xz = reshape(E_vec, Nx, Nz);

    % Extract x>=0 radial half line.
    mid = Nr;
    E_rz = E_xz(mid:end, :);
    I_rz = abs(E_rz).^2;

    info = struct();
    info.method = 'scalar 2D Helmholtz finite-difference frequency-domain';
    info.lambda = lambda;
    info.k0 = k0;
    info.x = x;
    info.A_size = size(A);
    info.nnz = nnz(A);
    info.assumptions = {'scalar TE/TM-like component','2D x-z model','not full vector Maxwell','boundaries affect result if too close'};
end
