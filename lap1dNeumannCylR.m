%% FINISHED

function L = lap1dNeumannCylR(r, dr)
    % 1D cylindrical laplacian calculator
    % ---------------------------------------------------------------------
    % Numerically calculates 1D laplacian (1/r*d(rd/dr))/dr) of r in 
    % cylindrical coordinatesb with Neumann BC (du/dr=0)
    % Calculation by:
    %         u"(i) ≈ (u(i−1)​−2u(i)​+u(i+1))/dr^2+(u(i+1)​​-u(i-1))/rdr​
    % =====================================================================
    % INPUTS:
    %        r - radial coordinate vector, Nr
    %        dr - distance between two grid points 
    % OUTPUT:
    %        L - laplacian approximation matrix 
    % *********************************************************************
    
    r = r(:);
    Nr = numel(r);

    main  = (-2/dr^2) * ones(Nr,1);
    upper = zeros(Nr-1,1);
    lower = zeros(Nr-1,1);

    % Interior nodes i = 2..Nr-1
    for i = 2:Nr-1
        lower(i-1) = 1/dr^2 - 1/(2*dr*r(i));
        upper(i)   = 1/dr^2 + 1/(2*dr*r(i));
    end
    
    
    main(1)  = -4/dr^2;
    upper(1) =  4/dr^2;

    main(Nr)     = -2/dr^2;
    lower(Nr-1)  =  2/dr^2;

    lower_full = [lower; 0];
    upper_full = [0; upper];

    L = spdiags([lower_full main upper_full], [-1 0 1], Nr, Nr);
end