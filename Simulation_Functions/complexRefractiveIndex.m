%% FIX absorption colorbar

function n = complexRefractiveIndex(pDiff, lambda)
    % Complex refractive index calculator
    % ---------------------------------------------------------------------
    % Evaluates the complex refractive index according to Bennett&Soref and
    % FCC diffusion for all space and time
    %                           nc = n + ik
    %                     k = alpha * lambda / (4*pi)
    %                            εr = nc^2
    % =====================================================================
    % INPUTS:
    %        pDiff - FCC distribution matrix, Nr x Nz x Nt, [1/m^3]
    %        lambda - beam's wavelength [m] 
    % OUTPUTS:
    %        n - complex refractive index matrix, Nr x Nz x Nt
    %        epsilon_r - relative permitivity matrix, Nr x Nz x Nt
    % *********************************************************************

    sp = systemParameters();
    
    [Nr,Nz,Nt] = size(pDiff);

    n = zeros(Nr, Nz, Nt);          % Complex refractive index, nc

    for it = 1:Nt
        N = pDiff(:,:,it);          
    
        % Bennett–Soref expects Ne, Nh at pump wavelength:
        [dn, dalpha] = BennettSoref(N, N); % Ne = Nh since 1 photon -> e + h
        
        n(:,:,it) = (sp.n + dn) + 1i*(sp.alpha + dalpha)*lambda/(4*pi);
    end
    
    % Calculating relative permitivity:
    %epsilon_r = n.^2;       % Assuming non-magnetic sample, μ = μ0
end