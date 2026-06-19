%% CHECKKKKKK

function [widthVal,fig] = PF_x(I, x, plotTitle, yTitle, xTitle)
    % Intensity plotter with automatic width marker
    % ---------------------------------------------------------------------
    % If profile is Gaussian-like:
    %   widthVal = regular FWHM
    %
    % If profile is donut-like:
    %   widthVal = distance between the two inner half-max points
    % ---------------------------------------------------------------------

    if nargin < 3 || isempty(plotTitle)
        plotTitle = '';
    end

   if nargin < 4 || isempty(yTitle)
        yTitle = 'Intensity [W/m^2]';
    end

    if nargin < 5 || isempty(xTitle)
        xTitle = 'x [\mum]';
    end

    % Force column vectors
    x = x(:);
    I = real(I(:));

    if length(x) ~= length(I)
        error("x and I must have the same length.");
    end

    % Remove small numerical negative values
    I(I < 0) = 0;

    x_um = x * 1e6;

    % Detect profile type
    isDonut = detectDonut(I, x);

    % Plot profile
    fig = figure('Color','w');
    plot(x_um, I, 'LineWidth', 2);
    hold on;
    axis tight; grid on;
    xlabel(xTitle, 'FontSize', 15);
    ylabel(yTitle, 'FontSize', 15);
    title(plotTitle);

    if isDonut
        widthVal = addDonutInnerDistance(I, x);
    else
        widthVal = addGaussianFWHM(I, x);
    end

    hold off;
end


%% ============================================================
%  Detect donut automatically
% ============================================================

function isDonut = detectDonut(I, x)

    leftIdx  = find(x < 0);
    rightIdx = find(x > 0);

    if isempty(leftIdx) || isempty(rightIdx)
        isDonut = false;
        return;
    end

    [~, centerIdx] = min(abs(x));

    leftMax  = max(I(leftIdx));
    rightMax = max(I(rightIdx));
    centerVal = I(centerIdx);

    sideMin = min(leftMax, rightMax);
    globalMax = max(I);

    % Donut condition:
    % 1. clear side lobes exist on both sides
    % 2. center is much lower than side lobes
    % 3. side lobes are significant compared to global max
    isDonut = centerVal < 0.5 * sideMin && ...
              sideMin > 0.4 * globalMax;
end


%% ============================================================
%  Gaussian FWHM
% ============================================================

function fwhmVal = addGaussianFWHM(I, x)

    Imax = max(I);
    halfMax = Imax / 2;

    if Imax <= 0
        warning("Cannot calculate FWHM: maximum intensity is zero.");
        fwhmVal = NaN;
        return;
    end

    [~, idxMax] = max(I);

    % Left half-max crossing
    idxLeft = find(I(1:idxMax) <= halfMax, 1, 'last');

    % Right half-max crossing
    idxRightLocal = find(I(idxMax:end) <= halfMax, 1, 'first');

    if isempty(idxLeft) || isempty(idxRightLocal)
        warning("Could not calculate Gaussian FWHM: profile does not cross half maximum on both sides.");
        fwhmVal = NaN;
        return;
    end

    idxRight = idxMax + idxRightLocal - 1;

    if idxLeft + 1 > length(I) || idxRight - 1 < 1
        warning("Invalid Gaussian FWHM indices.");
        fwhmVal = NaN;
        return;
    end

    % Interpolate crossings
    xHalfLeft = interpCrossing(x(idxLeft), x(idxLeft+1), ...
                               I(idxLeft), I(idxLeft+1), ...
                               halfMax);

    xHalfRight = interpCrossing(x(idxRight-1), x(idxRight), ...
                                I(idxRight-1), I(idxRight), ...
                                halfMax);

    fwhmVal = xHalfRight - xHalfLeft;

    % Plot FWHM line
    plot([xHalfLeft xHalfRight]*1e6, [halfMax halfMax], ...
         'LineWidth', 1.8);

    plot([xHalfLeft xHalfLeft]*1e6, [0 halfMax], ...
         ':', 'LineWidth', 1.2, color='r');

    plot([xHalfRight xHalfRight]*1e6, [0 halfMax], ...
         ':', 'LineWidth', 1.2, color='r');

    % plot(xHalfLeft*1e6, halfMax, 'o', ...
    %      'MarkerSize', 7, 'LineWidth', 2, color='r');
    % 
    % plot(xHalfRight*1e6, halfMax, 'o', ...
    %      'MarkerSize', 7, 'LineWidth', 2, color='r');

    text(mean([xHalfLeft xHalfRight])*1e6, halfMax, ...
        sprintf('%.3f \\mum', fwhmVal*1e6), ...
        'VerticalAlignment', 'bottom', 'FontSize',14);
end


%% ============================================================
%  Donut inner half-max distance
% ============================================================

function donutVal = addDonutInnerDistance(I, x)

    leftIdx  = find(x < 0);
    rightIdx = find(x > 0);

    if isempty(leftIdx) || isempty(rightIdx)
        warning("Cannot calculate donut distance: x must include negative and positive values.");
        donutVal = NaN;
        return;
    end

    % Find lobe peaks
    [IleftMax, localLeft] = max(I(leftIdx));
    [IrightMax, localRight] = max(I(rightIdx));

    idxLeftPeak  = leftIdx(localLeft);
    idxRightPeak = rightIdx(localRight);

    halfLeft  = IleftMax / 2;
    halfRight = IrightMax / 2;

    centerIdx = find(x >= 0, 1, 'first');

    % -----------------------------
    % Left inner half-max crossing
    % -----------------------------
    searchLeft = idxLeftPeak:centerIdx;
    idxCandidates = searchLeft(I(searchLeft) <= halfLeft);

    if isempty(idxCandidates)
        warning("Could not find inner half-max crossing for left donut lobe.");
        donutVal = NaN;
        return;
    end

    idxLeft2 = idxCandidates(1);
    idxLeft1 = idxLeft2 - 1;

    if idxLeft1 < 1
        warning("Invalid left donut crossing indices.");
        donutVal = NaN;
        return;
    end

    xHalfLeftInner = interpCrossing(x(idxLeft1), x(idxLeft2), ...
                                    I(idxLeft1), I(idxLeft2), ...
                                    halfLeft);

    % ------------------------------
    % Right inner half-max crossing
    % ------------------------------
    searchRight = centerIdx:idxRightPeak;
    idxCandidates = searchRight(I(searchRight) <= halfRight);

    if isempty(idxCandidates)
        warning("Could not find inner half-max crossing for right donut lobe.");
        donutVal = NaN;
        return;
    end

    idxRight1 = idxCandidates(end);
    idxRight2 = idxRight1 + 1;

    if idxRight2 > length(I)
        warning("Invalid right donut crossing indices.");
        donutVal = NaN;
        return;
    end

    xHalfRightInner = interpCrossing(x(idxRight1), x(idxRight2), ...
                                     I(idxRight1), I(idxRight2), ...
                                     halfRight);

    donutVal = xHalfRightInner - xHalfLeftInner;

    yLine = min(halfLeft, halfRight);

    % Plot lobe peaks
    % plot(x(idxLeftPeak)*1e6, IleftMax, 'o', ...
    %      'MarkerSize', 8, 'LineWidth', 2, color='k');
    % 
    % plot(x(idxRightPeak)*1e6, IrightMax, 'o', ...
    %      'MarkerSize', 8, 'LineWidth', 2, color='k');

    % Plot inner half-max points
    % plot(xHalfLeftInner*1e6, halfLeft, 'x', ...
    %      'MarkerSize', 9, 'LineWidth', 2, color='r');
    % 
    % plot(xHalfRightInner*1e6, halfRight, 'x', ...
    %      'MarkerSize', 9, 'LineWidth', 2, color='r');

    % Draw distance line
    plot([xHalfLeftInner xHalfRightInner]*1e6, [yLine yLine], ...
         '-', 'LineWidth', 1.8, color='r');

    plot([xHalfLeftInner xHalfLeftInner]*1e6, [0 yLine], ...
         ':', 'LineWidth', 1.2, color='r');

    plot([xHalfRightInner xHalfRightInner]*1e6, [0 yLine], ...
         ':', 'LineWidth', 1.2, color='r');

    text(mean([xHalfLeftInner xHalfRightInner])*1e6, yLine, ...
        sprintf('%.3f \\mum', donutVal*1e6), ...
        'VerticalAlignment', 'bottom', 'FontSize',14);
end


%% ============================================================
%  Linear interpolation helper
% ============================================================

function xCross = interpCrossing(x1, x2, y1, y2, yTarget)

    if y2 == y1
        xCross = mean([x1 x2]);
    else
        xCross = x1 + (yTarget - y1) * (x2 - x1) / (y2 - y1);
    end
end