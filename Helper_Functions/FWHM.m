%% FINISHED

function fwhm = FWHM(I, r)
    % FWHM calculator function
    % ---------------------------------------------------------------------
    % Calculates the FWHM of a given cartesian Gaussian profile for each
    % z according to:
    %                       FWHM(z) = 2*rhalf
    % where rhalf is the radius at half the peak
    % =====================================================================
    % INPUTS:
    %        I - beam's intensity profile, Nr x Nz [W/m^2]
    %        r - radial coordinate vector [m]
    % OUTPUT:
    %        fwhm - fwhm vector, Nz
    % *********************************************************************

    [Nr,Nz] = size(I);

    r = r(:);

    fwhm = zeros(Nz,1);

    for iz = 1:Nz
        Icol = I(:, iz);
        hm = 0.5 * max(Icol);

        % Central-lobe mask
        mask = (Icol >= hm);

        % Find end index of the connected region starting at r=0
        k_end = find(~mask, 1, 'first') - 1;

        if k_end >= Nr
            continue
        end

        % Interpolate crossing between k_end (>=hm) and k_end+1 (<hm):
        r1 = r(k_end); r2 = r(k_end + 1);
        I1 = Icol(k_end); I2 = Icol(k_end + 1);

        % Avoid division by zero:
        if I2 == I1
            rhalf = r1;
        else
            rhalf = r1 + (hm - I1) * (r2 - r1) / (I2 - I1);
        end

        fwhm(iz) = 2 * rhalf;
    end
end

% 2nd implementation using findpeaks:
% function fwhm_main = FWHM(I, x)
%     % FWHM calculator function
%     % ---------------------------------------------------------------------
%     % Calculates the main lobe's FWHM of a given cylindrical Gaussian
%     % profile for each z according to:
%     %                       FWHM(z) = 2*rhalf
%     % where rhalf is the radius at half the peak due to radial symmetry
%     %
%     % Also finds the width of other lobes if exist
%     % =====================================================================
%     % INPUTS:
%     %        I - beam's intensity profile, Nx x Nz [W/m^2]
%     %        x - x coordinate vector [m]
%     % OUTPUTS:
%     %        fwhm_main - main lobe FWHM vector, Nz [m]
%     % *********************************************************************
% 
%     x = x(:);
%     [~, Nz] = size(I);
% 
%     fwhm_main  = zeros(Nz,1);
% 
%     for iz = 1:Nz
%         Icol = I(:,iz);
% 
%         % Find all peaks and their widths at half peak:
%         [pks, ~, w] = findpeaks(Icol, x,'WidthReference', 'halfheight');
% 
%         if isempty(pks)
%             continue
%         end
% 
%         % Main lobe is the strongest peak:
%         [~, iMain] = max(pks);
% 
%         fwhm_main(iz) = w(iMain);   % Radial symmetry
%     end
% end
