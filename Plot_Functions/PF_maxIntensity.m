%% FINISHED

function PF_maxIntensity(I, t)
    % Compares intensity Plots
    % ---------------------------------------------------------------------
    % Plots two colormaps probe intensity before and after pump with
    % specific delay to check the affect on the probe beam
    % Note that we plot the maximal intensity ANYWHERE foor each time point
    % not at specific depth 
    % =====================================================================
    % INPUTS:
    %        I - intensity distribution matrix, Nx x Nz x Nt, [W/m^2]
    %        t - time vector [s]
    % *********************************************************************
    
    maxI = squeeze(max(I, [], [1 2]));
    %maxI = maxI/min(maxI);      % Normalization to original maximum

    figure;
    plot(t*1e12, maxI);
    axis tight; grid on;
    xlabel('t [ps]'); ylabel('Maximal Intensity [W/m^2]');
    title('Maximal Intensity');

end