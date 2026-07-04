%% FINISHED

function [D,r0,r1,r2] = computeDiameterTheoretically(z, type, l, lambda, w0, z0, n, plotTitle)
    % Calculates diameter theoretically
    % ---------------------------------------------------------------------
    % Calculates and plots the diameter of a given intensity profile matrix
    % (Nr x Nz, or Nx x Nz) according to:
    %                           D = 2w
    % where w is the length in r-axis where the intensity falls as 1/e^2.
    % =====================================================================
    % INPUTS:
    %        z - z coordinate, propagation vector [m]
    %        type - string ("Gauss" or "Donut")
    %        l - LG polynomial index
    %        lambda - beam's wavelength [m]
    %        w0 - waist radius at z=z0 [m]
    %        z0 - waist location along z [m]
    %        n - medium refractive index
    %        plotTitle - title, string (optional)
    % OUTPUTS:
    %        D - Diameter vector, Nz
    %        r0 - radius of donut maximum intensity, vector Nz
    %        r1, r2 - outer donut radii, vectors Nz
    % *********************************************************************
            
    if isempty(l)
        l = 1;
    end
    
    if nargin < 8 || ~(isstring(plotTitle))
        if string(type) == "Gauss"
            plotTitle = "Theoretical Gaussian Diameter";
        elseif string(type) == "Donut"
            plotTitle = "Theoretical Donut Diameters";
        else
            error("\nInappropriate type. Please use 'Gauss' or 'Donut'");
        end
    end

    zR = pi*n*w0^2/lambda;

    % Beam radius w(z)
    w = w0 * sqrt(1 + ((z-z0)/zR).^2);

    if string(type) == "Gauss"
        D = 2*w;
        FWHM = sqrt(2*log(2)) * w;

        figure;
        plot(z*1e6, D*1e6, 'k','LineWidth', 2, 'DisplayName', 'Diameter');
        hold on;
        plot(z*1e6, FWHM*1e6,'r','LineWidth', 2, 'DisplayName', 'FWHM');
        hold off;
        grid on; set(gcf,'Color','w');
        xlabel('z [\mum]', 'FontSize', 15);
        ylabel('Diameter [\mum]', 'FontSize', 15);
        title(plotTitle); legend('Location','best');
        
        r1 = 0; r2 = 0; r0 = 0;
    
    elseif string(type) == "Donut"
        x = -(exp(-1-2/l)*l)/l;
        
        % Inner and Outer radii:
        r0 = sqrt(abs(l)/2) .* w;
        r1 = sqrt(-l*lambertw(0,x)/2).*w;
        r2 = sqrt(-l*lambertw(-1,x)/2).*w;

        figure;
        plot(z*1e6, r1*2e6, 'k','LineWidth', 2, 'DisplayName', 'Inner diameter'); 
        hold on;
        plot(z*1e6, r2*2e6, 'r','LineWidth', 2, 'DisplayName', 'Outer diameter'); 
        plot(z*1e6, r0*2e6, 'b', 'LineWidth', 2, 'DisplayName', 'D_0');
        hold off; grid on; set(gcf,'Color','w');
        xlabel('z [\mum]', 'FontSize', 15);
        ylabel('Diameter [\mum]', 'FontSize', 15);
        title(plotTitle); legend('Location','best');

        D = 0;
    
    else
        error("\nInappropriate type. Please use 'Gauss' or 'Donut'");
    end
end