%% FINISHED

function prf = GB(r, z, lambda, w0, z0, E0)
    % Gaussian beam
    % ---------------------------------------------------------------------
    % Calculates the electric field (assuming polarizationi x direction),
    % or the spacial profile of a Gaussian beam
    % =====================================================================
    % INPUTS:
    %        r - radial coordinate vector [m]
    %        z - propagation vector [m]
    %        lambda - beam's wavelength [m]
    %        w0 - waist radius at z=z0 [m]
    %        E0 - amplitude at the origin, E(0,0)
    %        z0 - waist location along z [m]
    % OUTPUT:
    %        prf - complex field spacial profile in cylindrical coordinates
    % *********************************************************************
    
    Nr = length(r);
    Nz = length(z);

    prf = complex(zeros(Nr, Nz));

    % Beam constants:
    k  = 2*pi/lambda;
    zR = pi*w0^2/lambda;
    
    % Calculating E(x,y,z) for each element in z:
    for iz = 1:Nz
        % Shift the axis relative to the waist:
        Z = z(iz) - z0;
        
        % Beam radius w(z)
        w = w0 * sqrt(1 + (Z/zR)^2);
        
        % Gouy:
        gouy = atan(Z/zR);
        
        % Radius of curvature:
        if Z == 0
            R = Inf;
        else
            R = Z * (1 + (zR/Z)^2);
        end

        % Amplitude:
        amp = E0 * (w0/w) .* exp(-(r.^2)/(w^2));

        % Phase:
        phase = exp(-1i * (k*Z - gouy + k*(r.^2)/(2*R)));
        
        prf(:,iz) = amp .* phase;                
    end 
    
end