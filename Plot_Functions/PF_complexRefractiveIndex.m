%% FINISHED

function PF_complexRefractiveIndex(n_complex, r, z, x, lambda, iz, it, png)
    % Plots Different Refractive Index Plots
    % ---------------------------------------------------------------------
    % Plots colormaps of the imaginary and real parts of the complex
    % refractive index (n = n + ik) in xz plane for some t
    % 
    % Plots the refractive index and absorption coefficient [1/cm] at
    % certain z (iz) vs. x on the same graph
    % =====================================================================
    % INPUTS:
    %        n_complex - complex refractive index matrix at some t, Nr x Nz
    %        r - radial coordinate vector [m]
    %        z - z coordinate vector [m]
    %        x - x coordinate vector [m]
    %        lambda - beam's wavelength [m]
    %        iz - z index
    %        it - specific time point
    %        png - export to .png file (optional)
    % *********************************************************************
    
    arguments
        n_complex
        r
        z
        x
        lambda
        iz
        it
        png = false
    end

    sp = systemParameters();

    % Converting the matrix to cartesian:
    n_xz = cylToCart(n_complex, r, x);
    n_imag_x = 1e-2 * imag(n_xz) * 4 * pi/ lambda;  % α, [1/cm]
    
    % Refractive index:
    figure("Color",'w');
    imagesc(z*1e6, x*1e6, real(n_xz)); hold on
    xline(z(iz)*1e6, '--w', sprintf("z = %.2f \\mum", z(iz)*1e6), LineWidth=3);
    hold off
    
    set(gca,"YDir","normal"); axis tight;
    xlabel("z [\mum]", "FontSize",15);
    ylabel("x [\mum]", "FontSize",15);
    title(sprintf("Refraction at t=%d[ps]",it*1e12));
    colorbar;

    % Absorption coefficient:
    figure("Color",'w');
    imagesc(z*1e6, x*1e6, n_imag_x); hold on
    xline(z(iz)*1e6, '--w', sprintf("z = %.2f \\mum", z(iz)*1e6), LineWidth=3);
    hold off
    
    set(gca,"YDir","normal"); axis tight;
    xlabel("z [\mum]", "FontSize",15);
    ylabel("x [\mum]", "FontSize",15);
    title(sprintf("Absorption at t=%d[ps]",it*1e12));
    cb = colorbar(gca);
    cb.Label.String = "\alpha [1/cm]";
    
    % Combined graphs:
    hfig = figure("Color",'w');
    ax = gca;
    acolor = sp.randomColor();
    yyaxis right;  plot(x*1e6, n_imag_x(:,iz), 'Color', acolor, 'LineWidth',3, 'LineStyle',':');
    ylabel('$\alpha\ [\mathrm{cm}^{-1}]$', 'Interpreter','latex','FontSize',15);
    ax.YColor = acolor;%"#80B3FF"     
    
    yyaxis left; plot(x*1e6, real(n_xz(:,iz)), 'Color', [1 0.38 0.53], 'LineWidth', 3);
    ylabel('$n$','Interpreter', 'latex','FontSize', 15);
    ax.YColor = [1 0.38 0.53]; 

    axis tight; grid on; %box off;
    
    xlabel("x [$\mu$m]", "FontSize",15, 'Interpreter','latex');
    legend("n","$\alpha$", 'Interpreter','latex');
    title(sprintf('Refraction and absorption at $z=%.2f\\,\\mu\\mathrm{m}$ and $t=%.0f\\,\\mathrm{ps}$', ...
          z(iz)*1e6, it*1e12),'Interpreter', 'latex');    
    
    picturewidth = 25; hw_ratio = 0.65;
    
    set(findall(hfig,'-property','FontSize'), 'FontSize', 20);
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
    set(hfig, 'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);

    if png
        exportgraphics(hfig, "Figures_and_Results\Refraction and absorption.png", 'Resolution', 300);
    end

end