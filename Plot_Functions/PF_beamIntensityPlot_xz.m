%% FINISHED

function PF_beamIntensityPlot_xz(I, z, x, plotTitle)
    % Plots intensity in xz
    % ---------------------------------------------------------------------
    % Plots colormap of a beam's intensity in cartesian coordinates
    % =====================================================================
    % INPUTS:
    %        I - intensity spacial distribution matrix, Nx x Nz, [W/m^2]
    %        r - radial coordinate vector [m]
    %        z - z coordinate vector [m]
    %        x - x coordinate vector [m]
    %        plotTitle - part of the title, string
    % *********************************************************************
    
    if nargin < 4 || ~(isstring(plotTitle))
        plotTitle = '';
    end
    
    figure;
    imagesc(z*1e6, x*1e6, I);
    axis xy;
    xlabel('z [mum]');
    ylabel('x [mum]');
    title(plotTitle);
    colorbar;
    
end