%% FINISHED

function [dn,dalpha] = BennetSoref(Ne, Nh, lambda)
    % BennetSoref  Free-carrier refraction & absorption in Si
    % ---------------------------------------------------------------------
    % Implements empirical formulas (at 1550 nm):
    %               Δn = -8.8e-22*ΔNe - 8.5e-18*(ΔNh)^0.8
    %                   Δα =  8.5e-18*ΔNe + 6.0e-18*ΔNh
    % where ΔNe,ΔNh are in cm^-3 and Δα is in cm^-1.
    % =====================================================================
    % INPUT:
    %        Ne - electron concentration change [cm^-3]
    %        Nh - hole concentration change [cm^-3]
    %        lambda - wavelength [m]
    % OUTPUT:
    %        dn - Δn
    %        dalpha - Δα [1/m]
    % *********************************************************************

    % Empirical Soref-Bennett formulas at 1550 nm
    dn = -8.8e-22 .* Ne  - 8.5e-18 .* (Nh .^ 0.8);
    dalpha =  8.5e-16 .* Ne  + 6e-16 .* Nh;         % in [1/m]
    
    % Δk = Δα * λ / (4π), using λ in meters and α in 1/m
    dk = (dalpha .* lambda) / (4*pi);

end