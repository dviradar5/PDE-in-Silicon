%% CHECKKKKKKKKKK

function hfig = PF_probeComparison(Ilist, names, x, z, it, izMax, yTitle)
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
   
    if nargin < 7 || isempty(yTitle)
        yTitle = "Intensity [$\mathrm{W}/\mathrm{m}^2$]";
    end

    sp = systemParameters();

    % Number of intensities to compare:
    NI = numel(Ilist);

    figure('Color','w');
    for i = 1:NI
        subplot(NI,1,i);
        imagesc(z*1e6, x*1e6, Ilist{i});
        
        hold on;
        xline(z(izMax{i})*1e6, '--w', sprintf('z = %.3f \\mum', z(izMax{i})*1e6), LineWidth=3);
        hold off;

        axis xy; colorbar;
        
        ylabel('x [\mum]', 'FontSize',15);
        title(names{i});
    end
    xlabel('z [\mum]', 'FontSize',15);
    sgtitle(sprintf('Probe Intensity Maps at %.0f ps', it*1e12));

    % Comparing maximal intensity:
    hfig = figure('Color','w');
    hold on;

    styles = {'-', '--', ':', '-.'};
    %cmap = lines(NI);%cmap(i,:)

    for i = 1:NI
        I = Ilist{i};
        iz = izMax{i};

        style = styles{mod(i-1, numel(styles)) + 1};

        plot(x*1e6, I(:,iz), LineWidth=2, LineStyle=style, Color=sp.colors(i), ...
            DisplayName=sprintf('%s', names{i}));
        sprintf('%s  (z = %.2f \\mum)', names{i}, z(iz)*1e6);
    end

    hold off;
    
    axis tight; grid on;

    ax = gca;
    ax.FontSize = 18;

    xlabel("$x\,[\mu\mathrm{m}]$", "FontSize",15, 'Interpreter','latex');
    ylabel(yTitle, "FontSize",15, 'Interpreter','latex');
    legend('Location','best','Interpreter','latex');
    title('Maximal Probe Intensity Profiles');

    picturewidth = 25; hw_ratio = 0.65;
    
    set(findall(hfig,'-property','FontSize'), 'FontSize', 20);
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
    set(hfig, 'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);

end