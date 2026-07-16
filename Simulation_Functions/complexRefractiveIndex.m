%% FINISHED

function n = complexRefractiveIndex(pDiff, lambda)
    % Complex refractive index calculator
    % ---------------------------------------------------------------------
    % Evaluates the complex refractive index according to Bennett&Soref and
    % FCC diffusion for all space and time
    %                           nc = n + ik
    %                      k = α * lambda / (4*pi)
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
    
    [Nr,Nz] = size(pDiff);

    n = zeros(Nr, Nz);          % Complex refractive index, nc, Nt

    % Bennett–Soref expects Ne, Nh at pump wavelength of 1550nm:
    %[dn, dalpha] = BennettSoref(N, N); % Ne=Nh since 1 photon -> e + h

    % Fixing wavelength dependancy according to Drude model:
    %dn = dn .* (lambda / 1550e-9)^2;
    %dalpha = dalpha .* (lambda / 1550e-9)^2;

    % Drude model:
    [dn, dalpha] = Drude(pDiff, pDiff, lambda); % Ne=Nh since 1 photon -> e + h

    n = (sp.n + dn) + 1i*(sp.alpha + dalpha)*lambda/(4*pi);
    
    % Calculating relative permitivity, assuming non-magnetic sample μ=μ0:
    %epsilon_r = n.^2; 
end