%% FINISHED

function Ax = cylToCart(A, r, x)
    % Cylindrical to cartesian transformation
    % ---------------------------------------------------------------------
    % Transforms cylindrical vector A(:,z) to cartesian matrix, where A is
    % a symmetrical cylindrical matrix representing r and z
    % =====================================================================
    % INPUTS:
    %        A - matrix, Nr x Nz
    %        r - radial coordinate vector [m]
    %        x - x coordinate vector [m]
    % OUTPUT:
    %        Ax - spatial matrix in cartesian coordinates (xz), Nx x Nz
    % *********************************************************************
    
    r = r(:);
    x = x(:);

    [~,Nz] = size(A);
    
    R = abs(x);

    Ax = zeros(numel(x), Nz);

    for i = 1:Nz
        Ax(:,i) = interp1(r, A(:,i), R, 'linear', 0);
    end
end