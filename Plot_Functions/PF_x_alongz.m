%% FINISHED

function PF_x_alongz(I,x,z, plotTitle,ITitle,png)
    % Intensity plotter
    % ---------------------------------------------------------------------
    % Plots intensity vs. x at all z and specific t
    % =====================================================================
    % INPUTS:
    %        I - beam's intensity profile, Nx x Nz [W/m^2]
    %        x - x coordinate vector [m]
    %        z - z coordinate vector [m]
    %        plotTitle - part of the title, string (optional)
    %        ITitle - vertical axis' title, string (optional)
    %        png - export to .png file (optional)
    % *********************************************************************

    arguments
        I
        x
        z
        plotTitle = ""
        ITitle = "Intensity [W/m^2]"
        png = false
    end
    
    hfig = figure('Color','w');
    
    surf(z*1e6,x*1e6, I);
    shading interp; colorbar;
    
    axis tight; grid on; %box off;
    
    xlabel("z [$\mu$m]", 'FontSize',15, 'Interpreter','latex');
    ylabel("x [\mum]", 'FontSize',15, 'Interpreter','latex');
    zlabel(ITitle, 'FontSize',15, 'Interpreter','latex');
    title(plotTitle, 'Interpreter', 'latex');
    
    if png
        exportgraphics(hfig, "Figures_and_Results\" + plotTitle + ".png", 'Resolution', 300);
    end
end