%% FINISHED

function PF_intensity_x(I,x)
    % Intensity plotter
    % ---------------------------------------------------------------------
    % Plots intensity vs. x at specific z and t
    % =====================================================================
    % INPUTS:
    %        I - beam's intensity profile, Nx [W/m^2]
    %        x - x coordinate vector [m]
    % *********************************************************************

    figure;
    plot(x*1e6, I);
    axis tight; grid on;
    xlabel('x [\mum]'); ylabel('Intensity [W/m^2]');
    title('Intensity');
end