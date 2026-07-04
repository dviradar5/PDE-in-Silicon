%% FINISHED

function PF_x_alongz(I,x,z, plotTitle,ITitle)
    % Intensity plotter
    % ---------------------------------------------------------------------
    % Plots intensity vs. x at all z and specific t
    % =====================================================================
    % INPUTS:
    %        I - beam's intensity profile, Nx x Nz [W/m^2]
    %        x - x coordinate vector [m]
    %        z - z coordinate vector [m]
    %        plotTitle - part of the title, string
    %        ITitle - vertical axis' title, string
    % *********************************************************************

    if nargin < 4
        plotTitle = 'Intensity map along z';
    end
    
    if nargin < 5 || isempty(ITitle)
        ITitle = 'Intensity [W/m^2]';
    end
    
    figure('Color','w');
    surf(z*1e6,x*1e6, I);
    shading interp;
    colorbar;
    axis tight; grid on;
    xlabel('z [\mum]', 'FontSize',15);
    ylabel('x [\mum]', 'FontSize',15);
    zlabel(ITitle, 'FontSize',15);
    title(plotTitle, 'Interpreter', 'latex');
end