%% CHECKKKKKKKKKKKKKKKK

function h = interactiveIxMovie(Ixz, x, z, plotTitle, yLabel)
    % interactiveIxMovie
    % ---------------------------------------------------------------------
    % Creates an interactive "movie" of I(x) profiles along z.
    %
    % The input intensity matrix should be:
    %
    %       Ixz : Nx x Nz
    %
    % where each column Ixz(:,iz) is the transverse intensity profile
    % at z(iz).
    %
    % CONTROLS:
    %   Slider       - move between z indices
    %   Play/Pause   - automatic animation
    %   Prev/Next    - step manually
    %
    % Keyboard:
    %   Space        - play/pause
    %   Right arrow  - next z
    %   Left arrow   - previous z
    %
    % INPUTS:
    %   Ixz       - intensity matrix, Nx x Nz
    %   x         - transverse coordinate vector [m]
    %   z         - propagation coordinate vector [m]
    %   plotTitle - optional plot title
    %   yLabel    - optional y-axis label
    %
    % OUTPUT:
    %   h - structure of graphics handles
    % *********************************************************************

    if nargin < 4 || isempty(plotTitle)
        plotTitle = "Interactive I(x,z) Movie";
    end

    if nargin < 5 || isempty(yLabel)
        yLabel = "Intensity [W/cm^2]";
    end

    x = x(:);
    z = z(:).';

    [Nx, Nz] = size(Ixz);

    if numel(x) ~= Nx
        error("Length of x must match the number of rows in Ixz.");
    end

    if numel(z) ~= Nz
        error("Length of z must match the number of columns in Ixz.");
    end

    % Use real intensity only
    Ixz = real(Ixz);

    % Convert coordinates to microns for plotting
    x_um = x * 1e6;
    z_um = z * 1e6;

    % Initial index
    iz = 1;

    % Animation settings
    isPlaying = false;
    dt = 0.05;   % seconds between frames

    % ---------------------------------------------------------------------
    % Create figure and axes
    % ---------------------------------------------------------------------
    fig = figure( ...
        "Color", "w", ...
        "Name", plotTitle, ...
        "NumberTitle", "off", ...
        "KeyPressFcn", @keyPressCallback);

    ax = axes(fig, ...
        "Position", [0.10 0.25 0.85 0.68]);

    curve = plot(ax, x_um, Ixz(:,iz), "LineWidth", 1.8);
    
    levelLine = yline(ax, max(Ixz(:,iz))/exp(2), "r--", "I_{max}/e^2", ...
        "LineWidth", 1.5, ...
        "LabelHorizontalAlignment", "left", ...
        "HandleVisibility", "off");

    grid(ax, "on");
    xlabel(ax, "x [\mum]", "FontSize", 14);
    ylabel(ax, yLabel, "FontSize", 14);

    updateYLim(iz);

    title(ax, sprintf("%s | z = %.3f \\mum | index %d/%d", ...
        plotTitle, z_um(iz), iz, Nz), ...
        "FontSize", 14);

    % ---------------------------------------------------------------------
    % Slider step handling
    % ---------------------------------------------------------------------
    if Nz > 1
        sliderStep = [1/(Nz-1), min(10/(Nz-1), 1)];
    else
        sliderStep = [1, 1];
    end

    % ---------------------------------------------------------------------
    % UI controls
    % ---------------------------------------------------------------------
    slider = uicontrol(fig, ...
        "Style", "slider", ...
        "Units", "normalized", ...
        "Position", [0.10 0.12 0.65 0.05], ...
        "Min", 1, ...
        "Max", Nz, ...
        "Value", iz, ...
        "SliderStep", sliderStep, ...
        "Callback", @sliderCallback);

    playButton = uicontrol(fig, ...
        "Style", "pushbutton", ...
        "Units", "normalized", ...
        "Position", [0.78 0.12 0.08 0.05], ...
        "String", "Play", ...
        "Callback", @playPauseCallback);

    prevButton = uicontrol(fig, ...
        "Style", "pushbutton", ...
        "Units", "normalized", ...
        "Position", [0.10 0.05 0.08 0.05], ...
        "String", "Prev", ...
        "Callback", @prevCallback);

    nextButton = uicontrol(fig, ...
        "Style", "pushbutton", ...
        "Units", "normalized", ...
        "Position", [0.20 0.05 0.08 0.05], ...
        "String", "Next", ...
        "Callback", @nextCallback);

    indexText = uicontrol(fig, ...
        "Style", "text", ...
        "Units", "normalized", ...
        "Position", [0.32 0.05 0.35 0.05], ...
        "BackgroundColor", "w", ...
        "FontSize", 12, ...
        "String", sprintf("z = %.3f \\mum, index = %d/%d", ...
        z_um(iz), iz, Nz));

    % Return handles
    h.fig = fig;
    h.ax = ax;
    h.curve = curve;
    h.slider = slider;
    h.playButton = playButton;
    h.prevButton = prevButton;
    h.nextButton = nextButton;
    h.indexText = indexText;

    % =====================================================================
    % Nested callback functions
    % =====================================================================

    function updatePlot(newIndex)
        % Force index into legal range
        iz = round(newIndex);
        iz = max(1, min(Nz, iz));

        % Update curve
        set(curve, "YData", Ixz(:,iz));

        set(levelLine, "Value", max(Ixz(:,iz))/exp(2));

        % Update y-axis according to current profile
        updateYLim(iz);

        % Update slider
        set(slider, "Value", iz);

        % Update title
        title(ax, sprintf("%s | z = %.3f \\mum | index %d/%d", ...
            plotTitle, z_um(iz), iz, Nz), ...
            "FontSize", 14);

        % Update index text
        set(indexText, "String", ...
            sprintf("z = %.3f \\mum, index = %d/%d", ...
            z_um(iz), iz, Nz));

        drawnow;
    end

    function updateYLim(k)
        % Updates the y-axis limit according to the current frame only

        y = Ixz(:,k);
        ymax = max(y);

        if isempty(ymax) || ~isfinite(ymax) || ymax <= 0
            ymax = 1;
        end

        ylim(ax, [0, 1.05*ymax]);
    end

    function sliderCallback(src, ~)
        updatePlot(get(src, "Value"));
    end

    function playPauseCallback(~, ~)
        isPlaying = ~isPlaying;

        if isPlaying
            set(playButton, "String", "Pause");
        else
            set(playButton, "String", "Play");
            return;
        end

        while isPlaying && isvalid(fig)

            if iz < Nz
                updatePlot(iz + 1);
            else
                updatePlot(1);
            end

            pause(dt);
            drawnow;

            if ~isvalid(fig)
                break;
            end
        end
    end

    function prevCallback(~, ~)
        updatePlot(iz - 1);
    end

    function nextCallback(~, ~)
        updatePlot(iz + 1);
    end

    function keyPressCallback(~, event)
        switch event.Key
            case "space"
                playPauseCallback([], []);

            case "rightarrow"
                updatePlot(iz + 1);

            case "leftarrow"
                updatePlot(iz - 1);
        end
    end
end