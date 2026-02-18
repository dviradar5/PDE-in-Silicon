%% DOCUMENT

function fwhm = FWHM(I_prf, r)
    % FWHM calculator function
    % ---------------------------------------------------------------------
    % Calculates the FWHM of a given cylindrical Gaussian profile for each
    % z according to:
    %                       FWHM(z) = 2*rhalf
    % where rhalf is the radius at half the peak
    % =====================================================================
    % INPUTS:
    %        I_prf - beam's intensity profile, Nr x Nz [J]
    %        r - radial vector [m]
    % OUTPUT:
    %        fwhm - fwhm vector of laser_profile, Nz
    % *********************************************************************
       
    [Nr,Nz] = size(I_prf);
    
    r = r(:);

    fwhm = zeros(1,Nz);

    for iz = 1:Nz
        I = I_prf(:,iz);
        
        % Computing half of the peak:
        Imax = max(I);
        hm = 0.5*Imax;
        
        % Finding rhalf as the radius at which abs(I-hm) is minimal:
        [~, imin] = min(abs(I-hm));
        
        if (imin >= Nr)
            continue
        else
            fwhm(iz) = 2*r(imin);
        end
    end
end