%% FINISHED

function PF_colormapAnimation_xz(Arz, r, z, x, tIdx, vidName, cb_title, lockClim, fps)
    % Colormap animation
    % ---------------------------------------------------------------------
    % Plots and displays few frames one after the other in a finite loop.
    % Saves the video as .mp4 file.
    % =====================================================================
    % INPUTS:
    %        Arz - matrix, Nr x Nz
    %        r - radial spacial coordinate vector [m] 
    %        z - z coordinate vector [m] 
    %        x - x coordinate vector [m]
    %        tIdx - time index vector
    %        vidName - name of the video, string
    %        cb_title - title to put near the colorbar, string
    %        lockClim - boolean flag to lock the colorbar
    %        fps - number of frames per second
    % *********************************************************************

    zVec = z(:);
    xVec = x(:);
    
    if nargin < 8
        lockClim = true;
    end

    if nargin < 9
        fps = 1;
    end
    loops = 3;      % Repeat 3 times

    v = VideoWriter(vidName, 'MPEG-4');
    v.FrameRate = fps;
    open(v);

    fig = figure('Color','w');

    tl = tiledlayout(fig, 1, 1, 'Padding','compact', 'TileSpacing','compact');
    ax = nexttile(tl);

    % Converting to cartesian:
    Axz0 = cylToCart(Arz(:,:,tIdx(1)), r, x);

    hImg = imagesc(ax, zVec*1e6, xVec*1e6, Axz0);
    axis(ax, 'xy');
    axis(ax, 'tight');
    xlabel(ax, 'z [\mum]');
    ylabel(ax, 'x [\mum]');
    colormap(ax, 'parula');

    % Optionally lock colorbar limits:
    if lockClim
        climit = [min(Axz0(:)) max(Axz0(:))];
        if climit(1) == climit(2), climit = climit + [-1 1]; end
        clim(ax, climit);
    end
    

    % Same colorbar for all the frames:
    cb = colorbar(ax);
    cb.Label.String = cb_title;

    drawnow;

    for k = 1:loops
        for ii = 1:numel(tIdx)

            Axz = cylToCart(Arz(:,:,tIdx(ii)), r, x);
            set(hImg, 'CData', Axz);
            
            % Optionally unlock colorbar limits:
            if ~lockClim
                climit = [min(Axz(:)) max(Axz(:))];
                if climit(1) == climit(2), climit = climit + [-1 1]; end
                clim(ax, climit);
            end

            title(ax, sprintf('frame %d/%d', ii, numel(tIdx)));
            drawnow;

            writeVideo(v, getframe(fig));
            
            % Show one frame per second:
            pause(1/fps);
        end
    end

    close(v);
end
