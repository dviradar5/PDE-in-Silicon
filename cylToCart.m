%% FINISHED

function Axy = cylToCart(A, r, x)
    % Cylindrical to Cartesian transformation
    % ---------------------------------------------------------------------
    % Transforms cylindrical vector A(:,iz) to cartesian matrix. A is a
    % symmetrical cylindrical matrix representing r and z
    % =====================================================================
    % INPUTS:
    %        A - matrix, Nr x Nz
    %        r - radial coordinate vector [m]
    %        x - x coordinate vector [m]
    % OUTPUT:
    %        Pxy - spacial matrix in cartesian coordinates (xz), Nx x Nz
    % *********************************************************************
    
    r = r(:);
    x = x(:);

    [~,Nz] = size(A);
    
    R = abs(x);

    Axy = zeros(numel(x), Nz);

    for zz = 1:Nz
        Axy(:,zz) = interp1(r, A(:,zz), R, 'linear', 0);
    end
end