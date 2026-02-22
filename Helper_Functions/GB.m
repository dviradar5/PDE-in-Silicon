%% FINISHED

function prf = GB(r, z, lambda, w0, z0, E0)
    % Gaussian beam profile
    % ---------------------------------------------------------------------
    % Calculates the electric field (assuming polarization in x direction),
    % or the spatial profile of a Gaussian beam
    %
    % The expression of Gaussian beam is given by:
    %   E(r,z) = E0 * w0/w(z) * exp(-(r/w(z))^2 -i(k(z+r^2/2R(z))-gouy))
    % where gouy = arctan(z/zR)
    % =====================================================================
    % INPUTS:
    %        r - radial coordinate vector [m]
    %        z - z coordinate, propagation vector [m]
    %        lambda - beam's wavelength [m]
    %        w0 - waist radius at z=z0 [m]
    %        E0 - amplitude at the origin, E(0,0)
    %        z0 - waist location along z [m]
    % OUTPUT:
    %        prf - complex field spatial profile in cylindrical coordinates
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