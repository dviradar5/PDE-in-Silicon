function L = radialLaplacianMatrix(r)
%RADIALLAPLACIANMATRIX Axisymmetric transverse Laplacian matrix.
%
% Approximates
%       nabla_perp^2 E = d^2E/dr^2 + (1/r)dE/dr
% for a cylindrically symmetric field E(r), with regularity at r=0 and a
% zero-slope condition at r=rmax. A separate absorbing mask should be used
% near rmax for outgoing/edge suppression.

    r = r(:);
    Nr = numel(r);
    if Nr < 3
        error('Need at least 3 radial grid points.');
    end

    drs = diff(r);
    dr = mean(drs);
    if max(abs(drs - dr)) > 1e-9 * max(abs(dr), eps)
        error('This finite-difference radial Laplacian assumes uniform r spacing.');
    end

    L = spalloc(Nr, Nr, 3*Nr);

    if abs(r(1)) < 1e-12 * max(r(end), 1)
        % Regular axis: for even field, Laplacian(0) = 4*(E2-E1)/dr^2.
        L(1,1) = -4/dr^2;
        L(1,2) =  4/dr^2;
    else
        % If the grid does not include r=0, use the generic formula.
        ri = r(1);
        L(1,1) = -2/dr^2;
        L(1,2) =  2/dr^2; % weak one-sided substitute
        warning('r(1) is not zero; using approximate first-row boundary condition.');
    end

    for i = 2:Nr-1
        ri = r(i);
        L(i,i-1) = 1/dr^2 - 1/(2*dr*ri);
        L(i,i)   = -2/dr^2;
        L(i,i+1) = 1/dr^2 + 1/(2*dr*ri);
    end

    % Outer zero-slope boundary: dE/dr = 0 => ghost point E_{N+1}=E_{N-1}.
    L(Nr,Nr-1) =  2/dr^2;
    L(Nr,Nr)   = -2/dr^2;
end
