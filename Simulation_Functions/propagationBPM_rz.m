%% FINISHED

function [E_rz, I_rz] = propagationBPM_rz(E0_r, r, z, probe, n_complex, a, b)
    % Beam propagation calcualtor
    % ---------------------------------------------------------------------
    % Calaculates the probe's propagation inside the sample using BPM in
    % cylindrical coordinates with perfect matched layer boundary
    % conditions
    % =====================================================================
    % INPUTS:
    %        E0_r - complex electric field vector at z=0, Nr [V/m]
    %        r - radial coordinate vector, Nr [m]
    %        z - z coordinate, propagation vector [m]
    %        probe - probe laser beam, Laser-type object
    %        n_complex - complex refractive index matrix at specific time, Nr x Nz
    %        a - PML mask parameter
    %        b - PML mask parameter
    % OUTPUTS:
    %        E_rz - complex field inside the sample, Nr x Nz [V/m]
    %        I_rz - intensity in the sample, Nr x Nz [W/m^2]
    % *********************************************************************

    sp = systemParameters();

    rVec = r(:);

    [Nr, Nz] = size(n_complex);

    % Steps:
    dr = r(2)-r(1);
    dz = z(2)-z(1);

    k0 = 2*pi/probe.lambda;

    s = 1i * dz / (4*sp.n*k0);

    L = lap1dNeumannCylR(rVec, dr);
    Aimp = speye(Nr) - s*L;
    Aexp = speye(Nr) + s*L;
    
    % PML mask:
    mask = exp(-a * (rVec / rVec(end)).^b);
    
    E_rz = zeros(Nr, Nz);
    I_rz = zeros(Nr, Nz);
    
    E_rz(:,1) = E0_r(:);
    I_rz(:,1) = sp.eps0 * sp.c0 * real(n_complex(:,1)) .* (abs(E0_r(:)).^2)/2;

    E = E0_r(:);
    
    for iz = 2:Nz
        % Half-diffraction CN:
        E = Aimp \ (Aexp * E);

        % Medium step:
        n_real = real(n_complex(:,iz));
        n_imag = imag(n_complex(:,iz));

        % Phase from Δn:
        E = E .* exp( 1i * k0 * (n_real - sp.n) * dz );

        % Absorption from extinction coefficient:
        E = E .* exp( -k0 * n_imag * dz );

        % Half-diffraction CN:
        E = Aimp \ (Aexp * E);
        
        % PML mask to avoid boundary reflections:
        E = E .* mask;

        E_rz(:,iz) = E;
        
        I_rz(:,iz) = sp.eps0 * sp.c0 * n_real .* (abs(E).^2)/2;
    
    end
end
