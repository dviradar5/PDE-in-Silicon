%% FINISHED

function L = lap1dNeumann(N, d)
    % 1D cartesian laplacian calculator
    % ---------------------------------------------------------------------
    % Numerically calculates 1D laplacian (d^2/dx^2) of some cartesian
    % coordinate (say x) with Neumann BC (du/dx=0)
    % 
    % Calculation by:
    %                  u"(i) ≈(u(i−1)​−2u(i)​+u(i+1))/d^2​
    % =====================================================================
    % INPUTS:
    %        N - number of grid points
    %        d - distance between two grid points 
    % OUTPUT:
    %        L - laplacian approximation matrix 
    % *********************************************************************

    e = ones(N,1);
    S = spdiags([e -2*e e], [-1 0 1], N, N);    % Tridiagonal matrix

    % Neumann BC:
    S(1,2) = 2;
    S(N,N-1) = 2;

    L = S / d^2;
end
