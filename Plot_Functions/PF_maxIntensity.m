%% FINISHED

function PF_maxIntensity(I, t)
    % Calculates and plots maximal intensity
    % ---------------------------------------------------------------------
    % Plots the maximal intensity vs. time
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