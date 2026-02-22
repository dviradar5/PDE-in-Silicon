%% ADD LG EXPRESSION

function prf = LGB01(r, phi, z, lambda, w0, z0, E0)
    % Laguerre–Gaussian with p=0, l=1 (mode01)
    % ---------------------------------------------------------------------
    % Calculates the electric field (assuming polarization in x direction),
    % or the spatial profile of a Laguerre-Gaussian beam of mode 01
    % 
    % The expression of Laguerre-Gaussian beam of mode pl is given by:
    %  E(r,z) = E0 * w0/w(z) * (r*sqrt(2)/w(z))^|l| * Lp|l|(2(r/w(z))^2)
    %           * exp(-i(k(z+r^2/2R(z))+l*phi-(2p+|l|+1)gouy))
    % where Lp|l| is the generalized Lagguerre polynomial (including the
    % normalization constant) and gouy = arctan(z/zR)
    % =====================================================================
    % INPUTS:
    %        r - radial coordinate vector [m]
    %        phi - azimutal coordinate vector [rad]
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

    % Mode 01:
    l = 1;
    p = 0;

    % Beam constants:
    k  = 2*pi/lambda;
    zR = pi*w0^2/lambda;
    
    % Generalized Laguerre polynom 01:
    L01 = 1;
    
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
        amp = E0 * (w0/w) .* (sqrt(2)*r/w).^abs(l) .* L01 .* exp(-(r.^2)/(w^2));

        % Phase:
        phase = exp(1i * (k*Z - (abs(l)+1)*gouy + k*(r.^2)/(2*R)));
        vort_ph = exp(1i * l * phi);
        
        prf(:,iz) = amp .* vort_ph .* phase;                
    end
    
end