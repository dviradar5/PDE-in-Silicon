%% ADD FWHM LINE

function PF_x(I,x, plotTitle, verTitle, horTitle)
    % Intensity plotter
    % ---------------------------------------------------------------------
    % Plots intensity vs. x at specific z and t
    % =====================================================================
    % INPUTS:
    %        I - beam's intensity profile, Nx [W/m^2]
    %        x - x coordinate vector (or any other coordinate) [m]
    %        plotTitle - part of the title, string
    % *********************************************************************

    if nargin < 4
        verTitle = 'Intensity [W/m^2]';
        horTitle = 'x [\mum]';
    end
    
    figure;
    plot(x*1e6, I);
    axis tight; grid on;
    xlabel(verTitle); ylabel(horTitle);
    title(plotTitle);
end