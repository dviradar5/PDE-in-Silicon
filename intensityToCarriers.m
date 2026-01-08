%% FINISHED

function Ne = intensityToCarriers(I, alpha, lambda)
    % Convert intensity to carrier concentration
    % ---------------------------------------------------------------------
    % Calculates electron/hole concentration given intensity profile.
    % Assuming quantum efficiency of 1.
    % =====================================================================
    % INPUT:
    %        I - intensity spacial profile, Nr x Nz
    %        alpha - absorption coefficient [1/m] 
    %        lambda - beam's wavelength [m]
    % 
    % OUTPUT:
    %        Ne - carrier concentration [cm^-3]
    % *********************************************************************

    sp = systemParameters();
    
    % Single photon energy: 
    Eph = sp.h*sp.c0/lambda;        % hc/λ [J]
    
    F = I;          % fluence [J/m^2], since I isn't t-dependent (CHANGE)
    
    Ne  = alpha .* F / (Eph * 1e6); % [cm^-3]
end