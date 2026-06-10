%% FINISHED

function D = computeDiameter(I, r)
    % Calculates diameter
    % ---------------------------------------------------------------------
    % Calculates the diameter of a given intensity profile matrix (Nr x Nz,
    % or Nx x Nz) according to:
    %                           D = 2w
    % where w is the length in r-axis where the intensity falls as 1/e^2.
    % =====================================================================
    % INPUTS:
    %        I - intensity matrix, Nr x Nz
    %        r - coordinate vector
    % OUTPUT:
    %        D - Diameter vector, Nz
    % *********************************************************************

    [~,Nz] = size(I);
    
    for iz = 1:Nz
        prof = I(:,iz);

        [~, imax] = max(prof);

        prof_right = prof(imax:end);
        r_right = r(imax:end);

        idx = find(prof_right <= max(prof)/exp(2), 1, 'first');

        % Two points around crossing
        I1 = prof_right(idx-1); I2 = prof_right(idx);
        r1 = r_right(idx-1); r2 = r_right(idx);

        % Linear interpolation to find exact r where I = max/e^2
        w = r1 + (max(prof)/exp(2) - I1)*(r2 - r1)/(I2 - I1);

        D(iz) = 2*w;
    end
end