%% Document

function PF_combinePlots(L, labels, figs, yTitle, plotTitle)
    
    % New combined figure
    fAll = figure("Color",'w');
    axAll = axes(fAll);
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
                'LineWidth', 3, ...
                'DisplayName', string(oldName) + " - " + labels(i));
        end
    end
    
    hold(axAll, 'off');
    
    grid(axAll, 'on');
    xlabel(axAll, 'z [\mum]', 'FontSize', 25);
    ylabel(axAll, yTitle, 'FontSize', 25);
    %title(axAll, plotTitle);
    legend(axAll, 'Location', 'best');
end