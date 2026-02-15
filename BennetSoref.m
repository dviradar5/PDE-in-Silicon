%% FINISHED

function [dn,dalpha] = BennetSoref(Ne, Nh, lambda)
    % BennetSoref Free-carrier refraction & absorption in undoped Si
    % ---------------------------------------------------------------------
    % Implements empirical formulas (at 1550 nm):
    %               Δn = -8.8e-22*ΔNe - 8.5e-18*(ΔNh)^0.8
    %                   Δα =  8.5e-18*ΔNe + 6.0e-18*ΔNh
    % where ΔNe,ΔNh are in cm^-3 and Δα is in cm^-1.
    % =====================================================================
    % INPUT:
    %        Ne - electron concentration [cm^-3]
    %        Nh - hole concentration [cm^-3]
    %        lambda - wavelength [m]
    % OUTPUT:
    %        dn - Δn
    %        dalpha - Δα [1/m]
    % *********************************************************************
    
    sp = systemParameters();

    Ne0 = sp.ni;
    Nh0 = sp.ni;
    
    dNe = Ne - Ne0;
    dNh = Nh - Nh0;

    % Empirical Bennet-Soref formulas at 1550 nm
    dn = -8.8e-22 .* dNe  - 8.5e-18 .* (dNh .^ 0.8);
    dalpha =  8.5e-16 .* dNe  + 6e-16 .* dNh;         % in [1/m]
    
    % Δk = Δα * λ / (4π), using λ in meters and α in 1/m
    dk = (dalpha .* lambda) / (4*pi);

end