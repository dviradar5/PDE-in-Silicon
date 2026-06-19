%% FINISHED

function PF_beamIntensityPlot_xy(I, r, x, plotTitle)
    % Plots intensity in xy
    % ---------------------------------------------------------------------
    % Plots colormap of a beam's intensity in xy at specific z and time
    % =====================================================================
    % INPUTS:
    %        I - intensity radial distribution vector at some z, Nr, [W/m^2]
    %        r - radial coordinate vector [m]
    %        x - x coordinate vector [m]
    %        plotTitle - part of the title, string
    % *********************************************************************
    
    if nargin < 4 || isempty(plotTitle)
        plotTitle = 'Intensity Cross-Section';
    end

    [X,Y] = meshgrid(x,x);
    R = hypot(X,Y);
    Ixy = interp1(r, I, R, 'linear', 0);
    
    oom = log10(max(I));        % Order of magnitude
    zscale = 0.1 * 10^oom;

    figure('Color','w');
    surf(x*1e6, x*1e6, Ixy);
    shading interp;
    axis image; set(gca,'YDir','normal');
    xlabel('x [\mum]', 'FontSize',15); ylabel('y [\mum]', 'FontSize',15);
    title(plotTitle, 'Interpreter', 'latex');
    colorbar;
    daspect([1 1 zscale]);

end