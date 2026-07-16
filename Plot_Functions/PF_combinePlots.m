%% Document

function hfig = PF_combinePlots(L, labels, figs, yTitle, plotTitle)
    
    % New combined figure
    hfig = figure("Color",'w');
    axAll = axes(hfig);
    hold(axAll, 'on');
    
    figColors = lines(numel(L));
    
    for i = 1:numel(L)
    
        % Get the real plotting axes, not legend/colorbar axes
        axOld = findobj(figs(i), 'Type', 'axes');
        axOld = axOld(~strcmp(get(axOld, 'Tag'), 'legend'));
    
        % Get all lines from old axes
        linesOld = findobj(axOld, 'Type', 'line');
    
        % Reverse order so MATLAB preserves plotting order
        linesOld = flipud(linesOld);
    
        for j = 1:numel(linesOld)
            newLine = copyobj(linesOld(j), axAll);
    
            oldName = get(linesOld(j), 'DisplayName');  
    
            set(newLine, 'Color', figColors(i,:), ...
                'LineWidth', 3);
                %,DisplayName', string(oldName) + " - " + labels(i)
        end
    end
    
    hold(axAll, 'off');
    
    grid(axAll, 'on');
    xlabel(axAll, "$z\,[\mu\mathrm{m}]$", "FontSize",15, 'Interpreter','latex');
    ylabel(axAll, yTitle, "FontSize",15, 'Interpreter','latex');
    %title(axAll, plotTitle);
    legend(axAll, 'Location', 'best','Interpreter','latex');

    picturewidth = 25; hw_ratio = 0.65;
    
    set(findall(hfig,'-property','FontSize'), 'FontSize', 20);
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
    set(hfig, 'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);
end