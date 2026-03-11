%% FINISHED

function PF_plot_xz(I, z, x, plotTitle,izMax)
    % Plots matrix in xz
    % ---------------------------------------------------------------------
    % Plots colormap of a matrix in cartesian coordinates
    % =====================================================================
    % INPUTS:
    %        I - matrix, Nx x Nz, [W/m^2]
    %        z - z coordinate vector [m]
    %        x - x coordinate vector [m]
    %        plotTitle - part of the title, string
    %        izMax - z index at maximum intensity
    % *********************************************************************
    
    if nargin < 4 || ~(isstring(plotTitle))
        plotTitle = '';
    end
    
    figure;
    imagesc(z*1e6, x*1e6, I);
    
    if nargin == 5
        hold on
        xline(z(izMax)*1e6, '--w', sprintf('z = %.3f \\mum', z(izMax)*1e6), LineWidth=3);
        hold off
    end

    axis xy;
    xlabel('z [\mum]');
    ylabel('x [\mum]');
    title(plotTitle);
    colorbar;
    
end