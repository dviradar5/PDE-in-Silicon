%% DOCUMENT

function PF_complexRefractiveIndex(n_complex, r, z, x, lambda)
    % Plots:
    iz = 1;
    n_real_r = real(n_complex(:,iz,1));
    n_imag_r = imag(n_complex(:,iz,1));
    n_real_x = interp1(r, n_real_r, abs(x), 'linear', 'extrap');
    n_imag_x = interp1(r, n_imag_r, abs(x), 'linear', 'extrap');    % k
    n_imag_x = 1e-2 * n_imag_x * 4 * pi/ lambda;               % α, [1/cm]
    
    %%FIX AXIS%%%%%%%%%
    figure;
    imagesc(z*1e6, x*1e6, n_real_x);
    set(gca,'YDir','normal'); axis tight;
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title('reflection, t = 0 ns');
    colorbar;
    
    figure;
    imagesc(z*1e6, x*1e6, n_imag_x);
    set(gca,'YDir','normal'); axis tight;
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title('absorption, t = 0 ns');
    colorbar;
    
    % Combined graphs:
    figure;
    ax = gca;
    yyaxis left;  plot(x*1e6, n_imag_x, "r", LineWidth=3, LineStyle= ":");
    ylabel('\alpha [1/cm]'); ax.YColor = 'r';     
    yyaxis right; plot(x*1e6, n_real_x, "m", LineWidth=3);
    ylabel('n'); ax.YColor = 'm';
    axis tight; grid on;
    xlabel('x [\mum]');
    legend('imag(n)', 'real(n)');
    title('z=0m, t = 0 ns');

end