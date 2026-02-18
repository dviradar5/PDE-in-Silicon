%% DOCUMENT
%% SIMPLIFY
%% CHANGE NAME

function PF_FCCPlot_xz_from_rzt(pDiff, r, z, t, tIdx, fps)
%FCCAnimate_xz_loop  Animate pDiff(r,z,t) as p(x,z) over time, looping + slow.
%
%   FCCAnimate_xz_loop(pDiff, r, z)
%       20 frames, fps=2 (slow), loops forever on screen, no file saved.
%
%   FCCAnimate_xz_loop(pDiff, r, z, tIdx, gifPath, fps, nLoops, useInterp)
%       tIdx     : time indices to show (vector)
%       gifPath  : "" (no gif) or "file.gif"
%       fps      : frames per second (smaller = slower). Default 2.
%       nLoops   : number of on-screen loops. Default inf.
%       useInterp: true/false. Default false (tries fast exact mapping first).

    if nargin < 5 || isempty(tIdx)
        Nt = size(pDiff,3);
        tIdx = round(linspace(1, 100, 30));
        tIdx = unique(tIdx, 'stable');
    end
    if nargin < 6 || isempty(fps), fps = 2; end 

    rVec = r(:);
    zVec = z(:);
    [Nr, Nz, Nt] = size(pDiff);

    if numel(rVec) ~= Nr, error('length(r) must equal size(pDiff,1)'); end
    if numel(zVec) ~= Nz, error('length(z) must equal size(pDiff,2)'); end

    tIdx = tIdx(:).';
    if any(tIdx < 1 | tIdx > Nt | mod(tIdx,1)~=0)
        error('tIdx must contain integer indices in [1..Nt].');
    end

    % Symmetric x-grid: length 2*Nr-1, x=0 at center
    xVec = [-flipud(rVec(2:end)); rVec];
    absx = abs(xVec);

    % Try fast exact index mapping
    haveIdxMap = false;
    tol = 1e-12 * max(1, max(rVec));
    [tf, idxMap] = ismembertol(absx, rVec, tol);
    haveIdxMap = all(tf);


    % Fixed color limits across frames
    pMin = +inf; pMax = -inf;
    for k = 1:numel(tIdx)
        frameData = pDiff(:,:,tIdx(k));
        pMin = min(pMin, min(frameData(:)));
        pMax = max(pMax, max(frameData(:)));
    end
    if pMin == pMax, pMax = pMin + 1; end

    % Figure init
    fig = figure;
    ax = axes(fig);

    % First frame
    IxZ = toCartesianXZ(pDiff(:,:,tIdx(1)));
    hImg = imagesc(ax, zVec*1e6, xVec*1e6, IxZ);
    axis(ax, 'xy');
    xlabel(ax, 'z'); ylabel(ax, 'x');
    colorbar(ax);
    clim(ax, [pMin pMax]);

    loopCounter = 0;
    while isvalid(fig)
        loopCounter = loopCounter + 1;

        for i = 1:numel(tIdx)
            tt = tIdx(i);

            IxZ = toCartesianXZ(pDiff(:,:,tt));
            set(hImg, 'CData', IxZ);
            title(ax, sprintf('t index = %d   (frame %d/%d)', t(tt)*1e12, i, numel(tIdx)));

            drawnow;

            pause(1/fps);  % controls speed (slower when fps is small)
        end
    end

    % ---- helper ----
    function IxZ_local = toCartesianXZ(pr)
        if haveIdxMap
            IxZ_local = pr(idxMap, :);
        else
            IxZ_local = zeros(numel(xVec), Nz);
            for zz = 1:Nz
                IxZ_local(:,zz) = interp1(rVec, pr(:,zz), absx, 'linear', 'extrap');
            end
        end
    end
end