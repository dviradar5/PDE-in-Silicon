%% FINISHED

function [D,r1,r2] = computeDiameter(I, r, z, type,l, lambda, w0, z0)
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
    %        z - z coordinate, propagation vector [m]
    %        type - string ("Gauss" or "Donut")
    %        l - LG polynomial index
    %        lambda - beam's wavelength [m]
    %        w0 - waist radius at z=z0 [m]
    %        z0 - waist location along z [m]
    % OUTPUT:
    %        D - Diameter vector, Nz
    % *********************************************************************
            
    if nargin < 4 || isempty(l)
        l = 0;
    end

    sp = systemParameters();

    [~,Nz] = size(I);
    
    if string(type) == "Gauss"
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
        
        r1=0; r2=0;
    
    elseif string(type) == "Donut"
        x = -(exp(-1-2/l)*l)/l;
        
        zR = pi*sp.n*w0^2/lambda;
        
        % Beam radius w(z)
        w = w0 * sqrt(1 + ((z-z0)/zR).^2);
        
        % Inner and Outer radii:
        r1 = sqrt(-l*lambertw(0,x)/2).*w;
        r2 = sqrt(-l*lambertw(-1,x)/2).*w;
        
        D=0;
    
    else
        error("\nInappropriate type. Please use 'Gauss' or 'Donut'");
    end
end