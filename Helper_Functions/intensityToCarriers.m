%% CHECK ASSUMPTION

function N = intensityToCarriers(I, tau, lambda)
    % Converts intensity to carrier concentration
    % ---------------------------------------------------------------------
    % Calculates electron/hole concentration given the intensity profile
    % 
    % Assumes quantum efficiency of 1 (also means that Ne = Nh)
    % 
    % Fluence calculation:
    %              I(r,z,t) = I(r,z) * exp(−4ln(2)*(t/τ)^2​)
    %           => F = ∫​I(t)dt = I(r,z) * sqrt(pi/ln(2)) * τ/2
    % =====================================================================
    % INPUTS:
    %        I - intensity spatial profile, Nr x Nz
    %        tau - gaussian envelope duration [s] 
    %        lambda - beam's wavelength [m]
    % OUTPUT:
    %        N - carrier concentration [m^-3]
    % *********************************************************************

    sp = systemParameters();
    
    F = I * (tau/2) *sqrt(pi/log(2));       % [J/m^2]

    Eph = sp.h * sp.c0 / lambda;                % [J]

    N = sp.alpha .* F / Eph;                % [1/m^3]

end
