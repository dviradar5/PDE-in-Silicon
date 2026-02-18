%% FINISHED

function PF_probeComparison(Ibefore, Iafter, x, z, it)
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
    title(sprintf('Probe Intensity with Pump and %d[s] Delay',it));
    colorbar;

end