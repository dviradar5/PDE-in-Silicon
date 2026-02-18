%% FINISHED

function [dn,dalpha] = BennettSoref(Ne, Nh)
    % Bennett-Soref Free-carrier refraction & absorption in undoped Si
    % ---------------------------------------------------------------------
    % Implements Benett&Soref's empirical formulas (at 1550 nm):
    %               Δn = -8.8e-22*ΔNe - 8.5e-18*(ΔNh)^0.8
    %                  Δα =  8.5e-18*ΔNe + 6.0e-18*ΔNh
    % ΔNe, ΔNh are in 1/cm^3 and Δα is in 1/cm
    % =====================================================================
    % INPUTS:
    %        Ne - electron concentration [m^-3]
    %        Nh - hole concentration [m^-3]
    % OUTPUTS:
    %        dn - change in the real part of the refractive index
    %        dalpha - change in the imaginary part, [1/m]
    % *********************************************************************
    
    sp = systemParameters();

    Ne0 = sp.ni;
    Nh0 = sp.ni;
    
    dNe = Ne*1e-6 - Ne0;            % Converting Ne to 1/cm^3
    dNh = Nh*1e-6 - Nh0;            % Converting Nh to 1/cm^3

    % Empirical Bennett-Soref formulas at 1550 nm
    dn = -8.8e-22 .* dNe  - 8.5e-18 .* (dNh .^ 0.8);
    dalpha =  8.5e-16 .* dNe  + 6e-16 .* dNh;       % [1/m]
    
    % Δk calculation taking α in 1/m
    % dk = (dalpha .* lambda) / (4*pi);

end