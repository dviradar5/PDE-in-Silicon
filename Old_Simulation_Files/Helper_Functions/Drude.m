%% FINISHED

function [dn,dalpha] = Drude(Ne, Nh, lambda)
    % Drude free-carrier refraction & absorption in undoped Si
    % ---------------------------------------------------------------------
    % Implements Drude model formulas for the change in complex refractive
    % index:
    %        Δn = -(e*lambda/pi*c)^2*(ΔNe/me + ΔNh/mh)/(8ε0*nr)
    %  Δα = (e*lambda/pi*c)^2*e*(ΔNe/(mue*me*^2)+ΔNh/(muh*mh*^2)/(4c*ε0*nr)
    %                           mu = e*τ/m*
    % Where mh and me are the effective masses, e is the elementary charge,
    % nr is the original, real refractive index and mue and muh are the
    % charges' mobilities.
    % 
    % ΔNe, ΔNh are in 1/m^3 and Δα is in 1/m
    % =====================================================================
    % INPUTS:
    %        Ne - electron concentration [m^-3]
    %        Nh - hole concentration [m^-3]
    % OUTPUTS:
    %        dn - change in the real part of the refractive index
    %        dalpha - change in the imaginary part, [1/m]
    % *********************************************************************
    
    sp = systemParameters();

    Ne0 = sp.ni;    % In [1/cm^3]
    Nh0 = sp.ni;
    
    dNe = Ne - Ne0*1e6;            % Converting Ne to [1/m^3]
    dNh = Nh - Nh0*1e6;

    % Coefficient:
    A = ((sp.e*lambda/(2*pi*sp.c0))^2)/(sp.eps0*sp.n);
    
    % Mobility definition:
    mue = sp.e*sp.tau_D/sp.me_eff;
    muh = sp.e*sp.tau_D/sp.mh_eff;

    % Drude model formulas:
    dn = -A*(dNe/sp.me_eff + dNh/sp.mh_eff)/2;
    dalpha = sp.e*A*(dNe/(mue*sp.me_eff^2) + dNh/(muh*sp.mh_eff^2))/sp.c0;  % [1/m]
    %dalpha = sp.e*A*(dNe/(sp.mu_e*sp.me_eff^2) + dNh/(sp.mu_h*sp.mh_eff^2))/sp.c0;  % [1/m]
end