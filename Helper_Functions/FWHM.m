%% CHECKKKK

function fwhmVec = FWHM(I, r, peakWindow, minCentralFrac)
    % FWHM calculator for the CENTRAL radial lobe
    % ---------------------------------------------------------------------
    % Calculates the FWHM of the central lobe of a cylindrically symmetric
    % probe intensity profile I(r,z).
    %
    % This version avoids artificial jumps caused by side-lobes becoming
    % higher than the on-axis lobe.
    %
    % INPUTS:
    %   I              - radial intensity profile, Nr x Nz [W/m^2]
    %   r              - radial coordinate vector, Nr x 1 [m]
    %   peakWindow     - optional: radius near r=0 where central peak is
    %                    searched [m]. Default: 5*dr.
    %   minCentralFrac - optional: reject central lobe if its peak is less
    %                    than this fraction of global max. Default: 0.02.
    %
    % OUTPUT:
    %   fwhmVec        - central-lobe FWHM versus z, 1 x Nz [m]
    % ---------------------------------------------------------------------

    r = r(:);

    % Allow vector input
    if isvector(I)
        I = I(:);
    end

    [Nr, Nz] = size(I);

    if length(r) ~= Nr
        error("Length of r must match number of rows in I.");
    end

    if r(1) < -eps
        error("This FWHM function expects a radial coordinate r >= 0, not Cartesian x.");
    end

    dr = mean(diff(r));

    if nargin < 3 || isempty(peakWindow)
        peakWindow = 5 * dr;
    end

    if nargin < 4 || isempty(minCentralFrac)
        minCentralFrac = 0.02;
    end

    fwhmVec = NaN(1, Nz);

    peakSearchIdx = find(r <= peakWindow);

    if isempty(peakSearchIdx)
        peakSearchIdx = 1;
    end

    for iz = 1:Nz

        profile = real(I(:, iz));
        profile(~isfinite(profile)) = 0;
        profile(profile < 0) = 0;

        globalMax = max(profile);

        if globalMax <= 0
            continue;
        end

        % -------------------------------------------------------------
        % Central peak only:
        % Search only very close to r=0, not over the full profile.
        % This prevents side-lobes from becoming the measured peak.
        % -------------------------------------------------------------
        [centralPeak, relPeakIdx] = max(profile(peakSearchIdx));
        iPeak = peakSearchIdx(relPeakIdx);

        if centralPeak <= 0
            continue;
        end

        % Optional rejection: central lobe is too weak compared to side-lobes
        if centralPeak < minCentralFrac * globalMax
            fwhmVec(iz) = NaN;
            continue;
        end

        halfMax = 0.5 * centralPeak;

        % Find first right-side crossing after the central peak
        idxRel = find(profile(iPeak:end) <= halfMax, 1, 'first');

        if isempty(idxRel)
            continue;
        end

        i2 = iPeak + idxRel - 1;
        i1 = i2 - 1;

        if i1 < 1 || i2 > Nr
            continue;
        end

        r1 = r(i1);
        r2 = r(i2);

        I1 = profile(i1);
        I2 = profile(i2);

        if abs(I2 - I1) < eps
            rHalf = 0.5 * (r1 + r2);
        else
            rHalf = r1 + (halfMax - I1) * (r2 - r1) / (I2 - I1);
        end

        % Cylindrical symmetry: full width = diameter
        fwhmVec(iz) = 2 * rHalf;
    end
end