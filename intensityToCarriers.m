%% DOCUMENT

function N = intensityToCarriers(I, tau, lambda)
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
    %        N - carrier concentration [cm^-3]
    % *********************************************************************

    sp = systemParameters();

    F = peakI_to_fluence(I, tau);
    Eph = sp.h*sp.c0/lambda;           % [J]
    N = sp.alpha .* F / Eph;            % [1/m^3]


    % n_p=this.E/this.eV/this.h_bar*this.lambda/this.c0;
end

function F = peakI_to_fluence(Ipeak, tau)
    temporalInt = (tau/2)*sqrt(pi/log(2));  % [s]
    F = Ipeak * temporalInt;                % [J/m^2]
end
