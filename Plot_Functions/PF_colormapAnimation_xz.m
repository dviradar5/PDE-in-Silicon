%% DOCUMENT

function PF_colormapAnimation_xz(Arz, r, z, x, tIdx, vidName, cb_title)

    zVec = z(:);
    xVec = x(:);

    fps   = 1;
    loops = 3;

    v = VideoWriter(vidName, 'MPEG-4');
    v.FrameRate = fps;
    open(v);

    fig = figure('Color','w');

    % Make an axes area that reserves space (colorbar won't get clipped)
    tl = tiledlayout(fig, 1, 1, 'Padding','compact', 'TileSpacing','compact');
    ax = nexttile(tl);

    Axz0 = cylToCart(Arz(:,:,tIdx(1)), r, x);

    hImg = imagesc(ax, zVec*1e6, xVec*1e6, Axz0);
    axis(ax, 'xy');
    axis(ax, 'tight');
    xlabel(ax, 'z [\mum]');
    ylabel(ax, 'x [\mum]');
    colormap(ax, 'parula');

    % Lock color limits
    climit = [min(Axz0(:)) max(Axz0(:))];
    if climit(1) == climit(2), climit = climit + [-1 1]; end
    clim(ax, climit);

    % Colorbar attached to *this* axes, placed on the right
    cb = colorbar(ax);
    cb.Label.String = cb_title;

    drawnow;

    for k = 1:loops
        for ii = 1:numel(tIdx)

            Axz = cylToCart(Arz(:,:,tIdx(ii)), r, x);
            set(hImg, 'CData', Axz);

            title(ax, sprintf('frame %d/%d', ii, numel(tIdx)));
            drawnow;

            writeVideo(v, getframe(fig));
            pause(1/fps);
        end
    end

    close(v);
end
