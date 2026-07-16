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
        ITitle = "Intensity [$\mathrm{W}/\mathrm{m}^2$]"
        png = false
    end
    
    hfig = figure('Color','w');
    
    surf(z*1e6,x*1e6, I);
    
    shading interp; colorbar;
    axis tight; grid on; %box off;
    
    ax = gca;
    ax.FontSize = 18;

    xlabel("$z\,[\mu\mathrm{m}]$", "FontSize",15, 'Interpreter','latex');
    ylabel("$x\,[\mu\mathrm{m}]$", "FontSize",15, 'Interpreter','latex');
    zlabel(ITitle, 'FontSize',15, 'Interpreter','latex');
    title(plotTitle, 'Interpreter', 'latex');

    picturewidth = 25; hw_ratio = 0.65;
    
    set(findall(hfig,'-property','FontSize'), 'FontSize', 20);
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
    set(hfig, 'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);

    if png
        exportgraphics(hfig, "Figures_and_Results\" + plotTitle + ".png", 'Resolution', 300);
    end
end