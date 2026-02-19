%% COMBINE WITH SIMILAR PLOTTING FUNCTIONS

function PF_FWHM_z(fwhm, z, it)
    % FWHM plotter
    % ---------------------------------------------------------------------
    % Plots FWHM vs. 
    % =====================================================================
    % INPUTS:
    %        I - beam's intensity profile, Nr x Nz [J]
    %        z - z coordinate vector [m]
    %        it - specific time point 
    % *********************************************************************

    figure;
    plot(z*1e6, fwhm*1e6);
    axis tight; grid on;
    xlabel('z [\mum]'); ylabel('FWHM [\mum]');
    title(sprintf('FWHM at %d[ps]', it*1e12));
end