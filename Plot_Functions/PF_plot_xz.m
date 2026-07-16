%% FINISHED

function hfig = PF_plot_xz(I, z, x, plotTitle,izMax, png)
    % Plots matrix in xz
    % ---------------------------------------------------------------------
    % Plots colormap of a matrix in cartesian coordinates
    % =====================================================================
    % INPUTS:
    %        I - matrix, Nx x Nz, [W/m^2]
    %        z - z coordinate vector [m]
    %        x - x coordinate vector [m]
    %        plotTitle - part of the title, string (optinal)
    %        izMax - z index at maximum intensity (optinal)
    %        png - export to .png file (optional)
    % *********************************************************************
    
    arguments
        I
        z
        x
        plotTitle = ""
        izMax = NaN
        png = false
    end

    
    hfig = figure("Color",'w');
    imagesc(z*1e6, x*1e6, I);
    
    if ~isnan(izMax)
        hold on
        xline(z(izMax)*1e6, 'w--', sprintf("z = %.3f \\mum", z(izMax)*1e6), LineWidth=3);
        hold off
    end
    
    ax = gca;
    ax.FontSize = 18;

    xlabel("$z\,[\mu\mathrm{m}]$", "FontSize",15, 'Interpreter','latex');
    ylabel("$x\,[\mu\mathrm{m}]$", "FontSize",15, 'Interpreter','latex');
    title(plotTitle);
    
    colorbar; grid off; axis tight; %box off;
        
    picturewidth = 25; hw_ratio = 0.65;
    
    set(findall(hfig,'-property','FontSize'), 'FontSize', 20);
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
    set(hfig, 'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);
    
    if png
        exportgraphics(hfig, "Figures_and_Results\" + plotTitle + ".png", 'Resolution', 300);
    end
end