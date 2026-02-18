%% FINISHED

function PF_beamIntensityPlot_xy(I, r, x, iz)
    % Plots intensity in xy
    % ---------------------------------------------------------------------
    % Plots colormap of a beam's intensity in xy at specific z (iz) and
    % time (it)
    % =====================================================================
    % INPUTS:
    %        I - intensity radial distribution vector at z(iz), Nr, [W/m^2]
    %        r - radial coordinate vector [m]
    %        x - x coordinate vector [m]
    %        iz - z coordinate
    % *********************************************************************

    [X,Y] = meshgrid(x,x);
    R = hypot(X,Y);
    Ixy = interp1(r, I, R, 'linear', 0);
    
    oom = log10(max(I));        % Order of magnitude
    zscale = 0.1 * 10^oom;

    figure;
    surf(x*1e6, x*1e6, Ixy);
    shading interp;
    axis image; set(gca,'YDir','normal');
    xlabel('x [\mum]'); ylabel('y [\mum]');
    title(sprintf('Intensity at z=%d[m]', iz));
    colorbar;
    daspect([1 1 zscale]);

end