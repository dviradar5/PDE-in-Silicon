%% FINISHED

function PF_complexRefractiveIndex(n_complex, r, z, x, lambda, iz, it)
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
    % *********************************************************************
    
    % Converting the matrix to cartesian:
    n_xz = cylToCart(n_complex, r, x);
    n_imag_x = 1e-2 * imag(n_xz) * 4 * pi/ lambda;  % α, [1/cm]
    
    % Refractive index:
    figure;
    imagesc(z*1e6, x*1e6, real(n_xz));
    set(gca,'YDir','normal'); axis tight;
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title(sprintf('Refraction at t=%d[ps]',it*1e12));
    colorbar;

    % Absorption coefficient:
    figure;
    imagesc(z*1e6, x*1e6, n_imag_x);
    set(gca,'YDir','normal'); axis tight;
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title(sprintf('Absorption at t=%d[ps]',it*1e12));
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
    legend('\alpha', 'n');
    title(sprintf('Refraction and Absorption at z=%d[m] at t=%d[ps]',z(iz),it*1e12));

end