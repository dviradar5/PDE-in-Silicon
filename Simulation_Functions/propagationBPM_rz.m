%% Maybe change the n in the intensity because it changes
%% e^iomegat

function [E_rz, I_rz] = propagationBPM_rz(E0_r, r, z, probe, n_complex)
    % Beam propagation calcualtor
    % ---------------------------------------------------------------------
    % Calaculates the probe's propagation inside the sample using BPM in
    % cylindrical coordinates
    % =====================================================================
    % INPUTS:
    %        E0_r - complex electric field vector at z=0, Nr [V/m]
    %        r - radial coordinate vector, Nr [m]
    %        z - z coordinate, propagation vector [m]
    %        probe - probe laser beam, Laser-type object
    %        n_complex - complex refractive index matrix at specific time, Nr x Nz
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
    k = sp.n * k0;

    a = 1i * dz / (4*k);

    L = lap1dNeumannCylR(rVec, dr);
    Aimp = speye(Nr) - a*L;
    Aexp = speye(Nr) + a*L;


    E_rz = zeros(Nr, Nz);
    E_rz(:,1) = E0_r(:);

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

        E_rz(:,iz) = E;
    end

    I_rz = 0.5 * sp.n * sp.eps0 * sp.c0 * abs(E_rz).^2;

end
