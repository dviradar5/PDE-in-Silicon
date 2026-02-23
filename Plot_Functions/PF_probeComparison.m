%% FINISHED

function PF_probeComparison(Ibefore, Iafter, x, z, it, izMax)
    % Compares intensities before and after pump
    % ---------------------------------------------------------------------
    % Plots two colormaps probe intensity before and after pump with
    % specific delay to check the affect on the probe beam
    % =====================================================================
    % INPUTS:
    %        Ibefore - intensity distribution matrix before pump, Nx x Nz, [W/m^2]
    %        Iafter - intensity distribution matrix some time after pump, Nx x Nz, [W/m^2]
    %        x - x coordinate vector [m]
    %        z - z coordinate vector [m]
    %        it - specific time point
    %        izMax - specific z index where maximal focusing occur
    % *********************************************************************

    figure;
    subplot(2,1,1);
    imagesc(z*1e6, x*1e6, Ibefore);
    axis xy;
    ylabel('x [\mum]');
    title('Probe Intensity without Pump');
    colorbar;
    
    subplot(2,1,2);
    imagesc(z*1e6, x*1e6, Iafter);
    axis xy;
    xlabel('z [\mum]');
    ylabel('x [\mum]');
    title(sprintf('Probe Intensity with Pump and %d[ps] Delay',it*1e12));
    colorbar;
    
    % Comparing maximal intensity:
    figure;
    hold on;
    plot(x*1e6, Ibefore(:,1), "r", LineWidth=2, LineStyle= ":");   
    plot(x*1e6, Iafter(:,izMax), "m", LineWidth=2);
    hold off;
    axis tight; grid on;
    xlabel('x [\mum]');
    ylabel('Intensity [W/m^2]');
    legend('before','after');
    title('Maximal Probe Intensity with and without Pump');

end