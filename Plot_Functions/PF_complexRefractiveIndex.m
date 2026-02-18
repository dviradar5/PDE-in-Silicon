%% FINISHED

function PF_complexRefractiveIndex(n_complex, r, z, x, lambda, iz)
    % Plots Different Refractive Index Plots
    % ---------------------------------------------------------------------
    % Plots colormaps of the imaginary and real parts of the complex
    % refractive index (n = n + ik) in xz plane
    % 
    % Plots the refractive index and absorption coefficient [1/cm] at
    % certain z (iz) vs. x on the same graph
    % =====================================================================
    % INPUTS:
    %        n_complex - complex refractive index matrix, Nr x Nz
    %        r - radial coordinate vector [m]
    %        z - z coordinate vector [m]
    %        x - x coordinate vector [m]
    %        lambda - beam's wavelength [m]
    %        iz - z index
    % *********************************************************************
    
    % Converting the matrix to cartesian:
    n_xz = cylToCart(n_complex, r, x);
    n_imag_x = 1e-2 * imag(n_xz) * 4 * pi/ lambda;  % α, [1/cm]
    
    % Refractive index:
    figure;
    imagesc(z*1e6, x*1e6, real(n_xz));
    set(gca,'YDir','normal'); axis tight;
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title('reflection');
    colorbar;

    % Absorption coefficient:
    figure;
    imagesc(z*1e6, x*1e6, n_imag_x);
    set(gca,'YDir','normal'); axis tight;
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title('absorption');
    cb = colorbar(gca);
    cb.Label.String = "\alpha [1/cm]";
    
    % Combined graphs:
    figure;
    ax = gca;
    yyaxis right;  plot(x*1e6, n_imag_x(:,iz), "r", LineWidth=3, LineStyle= ":");
    ylabel('\alpha [1/cm]'); ax.YColor = 'r';     
    yyaxis left; plot(x*1e6, real(n_xz(:,iz)), "m", LineWidth=3);
    ylabel('n'); ax.YColor = 'm';
    axis tight; grid on;
    xlabel('x [\mum]');
    legend('~imag(n)', 'real(n)');
    title(sprintf('Refraction and Absorption at z=%d[m]',z(iz)));

end