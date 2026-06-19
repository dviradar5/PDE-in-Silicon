%% FINISHED

function Ith = damageThreshold(I, tau, wl, title)
    % Calculates damage-threshold
    % ---------------------------------------------------------------------
    % Calculates the intensity damage-threshold for silicon illuminated by
    % a laser beam with specific pulse width and wavelength below 1[micron]
    %
    % Threshold flounce value was taken from Wang et al. (2010) (see paper
    % 033103_1_online in the papers folder)
    % =====================================================================
    % INPUTS:
    %        I - intensity matrix, Nr x Nz, [W/m^2]
    %        tau - laser pulse width, [s]
    %        wl - wavelength, [m]
    %        title - string
    % OUTPUT:
    %        Ith - damage threshold maximal intensity, [W/cm^2]
    % *********************************************************************
    
    % Damage-threshold flounce for laser with 775[nm] wavelength:
    Fth_10ps = 0.7;                                 % [J/cm^2]
    Fth = Fth_10ps*sqrt(tau/10e-12)*(wl/1064e-9);   % Corrections

    Ith = Fth/tau;                                  % [W/cm^2]
    
    [xMax,zMax,~] = findMax(I);
    Imax = I(xMax,zMax) * 1e-4;

    fprintf("\n%s", title);

    if Imax >= Ith
        fprintf("\nCaution!" + ...
            "\nyour intensity has passed the damage threshold " + ...
            "for Silicon with laser beam with wavelength of %.0f[nm] and " + ...
            "pulse width of %.0f[ps]", wl*1e9, tau*1e12);
    else
        fprintf("\nAll good ;)");
    end
    
    fprintf("\nMaximal intensity reached: %.3e[W/cm^2]", Imax);
    fprintf("\nDamage-threshold maximal intensity: %.3e[W/cm^2]\n", Ith);

end