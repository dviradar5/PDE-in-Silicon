%% FINISHED

function [Pxy, xPlot, yPlot] = cylToCart(A, r, iz)
    % Cylindrical to Cartesian transformation
    % ---------------------------------------------------------------------
    % Transforms cylindrical vector A(:,iz) to cartesian matrix. A is a
    % symmetrical cylindrical matrix representing r and z
    % =====================================================================
    % INPUTS:
    %        A - matrix, Nr x Nz
    %        r - radial coordinate vector [m]
    %        iz - z index (constant z) [m]
    % OUTPUT:
    %        xplot - x vector [m]
    %        yplot - y vector [m]
    %        prf - complex field spacial profile in cylindrical coordinates
    % *********************************************************************
    
    rmax = max(r);
    Nr = length(r);

    xPlot = linspace(-rmax, rmax, 2*Nr-1);
    yPlot = xPlot;                          % Radial symmetry
    
    r = r(:);
    
    p_r = A(:, iz);

    [X,Y] = meshgrid(xPlot, yPlot);
    R = hypot(X,Y);

    Pxy = interp1(r, p_r(:), R, 'linear', 0);   % Ny x Nx
end