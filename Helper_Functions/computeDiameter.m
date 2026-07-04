%% CHECKKKKKKKKKKKKK

function [fig,D,r0,r1,r2] = computeDiameter(I, r, z, type, plotTitle)
    % Calculates diameter
    % ---------------------------------------------------------------------
    % Calculates and plots the diameter of a given intensity profile matrix
    % (Nr x Nz, or Nx x Nz) using interpolation
    % =====================================================================
    % INPUTS:
    %        I - beam's intensity profile, Nr x Nz [W/m^2]
    %        r - radial coordinate vector [m]
    %        z - z coordinate, propagation vector [m]
    %        type - string ("Gauss" or "Donut")
    %        plotTitle - title, string (optional)
    % OUTPUTS:
    %        fig - figure object created
    %        D - Diameter vector, Nz
    %        r0 - radius of donut maximum intensity, vector Nz
    %        r1, r2 - outer donut radii, vectors Nz
    % *********************************************************************
    
    
    if nargin < 5 || ~(isstring(plotTitle))
        if string(type) == "Gauss"
            plotTitle = "Gaussian Diameter";
        elseif string(type) == "Donut"
            plotTitle = "Donut Diameters";
        else
            error("\nInappropriate type. Please use 'Gauss' or 'Donut'");
        end
    end

    r = r(:);
    z = z(:).';

    [~,Nz] = size(I);

    isSignedCoordinate = (min(r) < 0) && (max(r) > 0);

    if string(type) == "Gauss"
        D = nan(1,Nz);
        FWHM = nan(1,Nz);

        for iz = 1:Nz
            prof = real(I(:,iz));
            prof = prof(:);

            [Imax, imax] = max(prof);
            
            %D(iz) = interpWidthAtLevel(prof, r, imax, Imax/exp(2), isSignedCoordinate);
            FWHM(iz) = interpWidthAtLevel(prof, r, imax, Imax/2, isSignedCoordinate);
        end
            
        D = sqrt(2/log(2)) .* FWHM;

        fig = figure("Color",'w');
        plot(z*1e6, D*1e6, 'k','LineWidth', 2, 'DisplayName', 'Diameter');
        hold on;
        plot(z*1e6, FWHM*1e6,'k--','LineWidth', 2, 'DisplayName', 'FWHM');
        hold off; grid on;
        xlabel('z [\mum]', 'FontSize', 15);
        ylabel('Diameter [\mum]', 'FontSize', 15);
        title(plotTitle); legend('Location','best');

        r0 = 0; r1 = 0; r2 = 0;

    elseif string(type) == "Donut"
        r0 = nan(1,Nz);
        r1 = nan(1,Nz);
        r2 = nan(1,Nz);

        % Donut radius should be measured on a radial coordinate.
        % If a signed x coordinate is supplied, use only the x >= 0 side.
        if isSignedCoordinate
            keep = r >= 0;
            r_donut = r(keep);
            I_donut = I(keep,:);
        else
            r_donut = r;
            I_donut = I;
        end

        for iz = 1:Nz
            prof = real(I_donut(:,iz));
            prof = prof(:);

            [Imax, imax] = max(prof);

            if Imax <= 0 || isnan(Imax)
                continue;
            end

            target = Imax/exp(2);

            % Radius of maximum ring intensity
            r0(iz) = r_donut(imax);

            % Inner 1/e^2 crossing, before the maximum ring
            idxInner = find(prof(1:imax) <= target, 1, 'last');
            if ~isempty(idxInner) && idxInner < imax
                r1(iz) = interpCrossing( ...
                    r_donut(idxInner),   prof(idxInner), ...
                    r_donut(idxInner+1), prof(idxInner+1), ...
                    target);
            end

            % Outer 1/e^2 crossing, after the maximum ring
            idxOuterLocal = find(prof(imax:end) <= target, 1, 'first');
            if ~isempty(idxOuterLocal) && idxOuterLocal > 1
                idxOuter = imax + idxOuterLocal - 1;
                r2(iz) = interpCrossing( ...
                    r_donut(idxOuter-1), prof(idxOuter-1), ...
                    r_donut(idxOuter),   prof(idxOuter), ...
                    target);
            end
        end

        fig = figure("Color",'w');
        plot(z*1e6, r1*2e6, 'k','LineWidth', 2, 'DisplayName', 'Inner diameter');
        hold on;
        plot(z*1e6, r2*2e6, 'k--','LineWidth', 2, 'DisplayName', 'Outer diameter');
        plot(z*1e6, r0*2e6, 'k:', 'LineWidth', 2, 'DisplayName', 'D_0');
        hold off; grid on;
        xlabel('z [\mum]', 'FontSize', 15);
        ylabel('Diameter [\mum]', 'FontSize', 15);
        title(plotTitle); legend('Location','best');
        
        D = 0;

    else
        error("\nInappropriate type. Please use 'Gauss' or 'Donut'");
    end
end

% =====================================================================
% Helper Functions:
% =====================================================================

function width = interpWidthAtLevel(prof, coord, imax, target, isSignedCoordinate)
    % Returns full width at a chosen intensity level.
    %
    % If coord is signed Cartesian x:
    %   width = xRight - xLeft.
    %
    % If coord is radial r >= 0:
    %   width = 2*rCross.

    width = nan;

    if isSignedCoordinate
        % Left crossing
        idxLeft = find(prof(1:imax) <= target, 1, 'last');
        if isempty(idxLeft) || idxLeft >= imax
            return;
        end

        xLeft = interpCrossing( ...
            coord(idxLeft),   prof(idxLeft), ...
            coord(idxLeft+1), prof(idxLeft+1), ...
            target);

        % Right crossing
        idxRightLocal = find(prof(imax:end) <= target, 1, 'first');
        if isempty(idxRightLocal) || idxRightLocal <= 1
            return;
        end

        idxRight = imax + idxRightLocal - 1;

        xRight = interpCrossing( ...
            coord(idxRight-1), prof(idxRight-1), ...
            coord(idxRight),   prof(idxRight), ...
            target);

        width = xRight - xLeft;

    else
        % Radial coordinate: find the crossing to the right of the peak.
        idxRightLocal = find(prof(imax:end) <= target, 1, 'first');
        if isempty(idxRightLocal) || idxRightLocal <= 1
            return;
        end

        idxRight = imax + idxRightLocal - 1;

        rCross = interpCrossing( ...
            coord(idxRight-1), prof(idxRight-1), ...
            coord(idxRight),   prof(idxRight), ...
            target);

        width = 2*rCross;
    end
end

function xCross = interpCrossing(x1, y1, x2, y2, yTarget)
    % Linear interpolation for the coordinate where y = yTarget.

    if y2 == y1
        xCross = (x1 + x2)/2;
    else
        xCross = x1 + (yTarget - y1)*(x2 - x1)/(y2 - y1);
    end
end
