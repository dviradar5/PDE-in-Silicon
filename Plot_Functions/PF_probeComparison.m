%% FINISHED

function PF_probeComparison(Ilist, names, x, z, it, izMax)
    % Compares intensities
    % ---------------------------------------------------------------------
    % Plots several colormaps of (probe) intensity before and after pump
    % to check the effect on the probe beam
    % =====================================================================
    % INPUTS:
    %        Ilist - list of intensity distribution matrices, where the
    %                first one is before the pump, Nx x Nz, [W/m^2]
    %        names - list of intensity titles (for the legend)
    %        x - x coordinate vector [m]
    %        z - z coordinate vector [m]
    %        it - specific time point
    %        izMax - list of z index where maximal focusing occur
    % *********************************************************************
    
    % Number of intensities to compare:
    NI = numel(Ilist);

    figure;
    for i = 1:NI
        subplot(NI,1,i);
        imagesc(z*1e6, x*1e6, Ilist{i});
        
        hold on;
        xline(z(izMax{i})*1e6, '--w', sprintf('z = %.3f \\mum', z(izMax{i})*1e6), LineWidth=3);
        hold off;

        axis xy;
        colorbar;
        
        xlabel('z [\mum]');
        ylabel('x [\mum]');
        title(names{i});
    end
    sgtitle(sprintf('Probe Intensity Maps at %.0f ps', it*1e12));

    % Comparing maximal intensity:
    figure;
    hold on;

    styles = {'-', '--', ':', '-.'};
    cmap = lines(NI);

    for i = 1:NI
        I = Ilist{i};
        iz = izMax{i};

        style = styles{mod(i-1, numel(styles)) + 1};

        plot(x*1e6, I(:,iz), LineWidth=2, LineStyle=style, Color=cmap(i,:), ...
            DisplayName=sprintf('%s  (z = %.2f \\mum)', names{i}, z(iz)*1e6));
    end

    hold off;
    axis tight; grid on;
    xlabel('x [\mum]');
    ylabel('Intensity [W/m^2]');
    legend('Location','best');
    title('Maximal Probe Intensity Profiles');

end