%% FINISHED

function PF_FWHM_z(fwhm, z, plotTitle)
    % FWHM plotter
    % ---------------------------------------------------------------------
    % Plots FWHM vs. z. Also points out the minimum FWHM
    % =====================================================================
    % INPUTS:
    %        I - beam's intensity profile, Nr x Nz [J]
    %        z - z coordinate vector [m]
    %        plotTitle - plot title, string 
    % *********************************************************************
    
    if nargin < 3
        plotTitle = 'FWHM';
    end
    
    [minFWHM,i] = min(fwhm);
    fprintf('Minimum FWHM: %.4e m, Achieved at z = %.2f um\n', minFWHM, z(i)*1e6);

    figure('Color','w');
    plot(z*1e6, fwhm*1e6,HandleVisibility="off");
    hold on;
    plot(z(i)*1e6,minFWHM*1e6,'x',DisplayName="Minimum", MarkerSize=10);
    hold off;
    
    grid on;
    xlabel('z [\mum]', 'FontSize',15); ylabel('FWHM [\mum]', 'FontSize',15);
    title(plotTitle);
    legend show;
end