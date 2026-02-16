%% DOCUMENT

function PF_beamPlotXY(r, x, Ibeam, zi, ti)
    [X,Y] = meshgrid(x,x);
    R = hypot(X,Y);
    Ixy = interp1(r, Ibeam, R, 'linear', 0);

    figure;
    surf(x*1e6, x*1e6, Ixy);
    shading interp;
    axis image; set(gca,'YDir','normal');
    xlabel('x [\mum]'); ylabel('y [\mum]');
    title(sprintf('Intensity at z=%.3f mum, t=%.3f ns', zi*1e6, ti*1e9)); %%%%%%%%%%%%%%%%%%%writes ti instead of t(i)
    colorbar;
end