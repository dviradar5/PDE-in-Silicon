%% Maybe change the n in the intensity because it changes

function [E_rz, I_rz] = propagationBPM_rz(E0_r, r, z, probe, n_complex_rz)
    % Beam propagation calcualtor
    % ---------------------------------------------------------------------
    % Calaculates the probe's propagation inside the sample using BPM in
    % cylindrical coordinates
    % =====================================================================
    % INPUTS:
    %        E0_r - complex electric field vector at z=0 [V/m]
    %        r - radial coordinate vector, Nr [m]
    %        z - z coordinate, propagation vector [m]
    %        probe - probe laser beam, Laser-type object
    %        n_complex_rz - complex refractive index matrix at specific time, NrxNz
    % OUTPUTS:
    %        E_rz - complex field, NrxNz [V/m]
    %        I_rz - intensity, NrxNz [W/m^2]
    % *********************************************************************

    sp = systemParameters();

    rVec = r(:);

    [Nr, Nz] = size(n_complex_rz);

    % Steps:
    dr = r(2)-r(1);
    dz = z(2)-z(1);

    k0 = 2*pi/probe.lambda;
    k = sp.n * k0;

    % Build axisymmetric radial Laplacian operator L such that:
    %   (L E)_i ≈ d2E/dr2 + (1/r) dE/dr   (m=0)
    %
    % CN half-step for diffraction:
    %   (I - a L) E^{*} = (I + a L) E
    % where a = i dz/(4k)
    a = 1i * dz / (4*k);

    L = lap1dNeumannCylR(rVec, dr);
    Aimp = speye(Nr) - a*L;
    Aexp = speye(Nr) + a*L;


    E_rz = zeros(Nr, Nz);
    E_rz(:,1) = E0_r(:);

    E = E0_r(:);

    for iz = 2:Nz
        % ---- half diffraction (CN) ----
        E = Aimp \ (Aexp * E);

        % ---- medium step (phase + absorption) ----
        n_real = real(n_complex_rz(:,iz));
        n_imag = imag(n_complex_rz(:,iz));

        % phase from Δn relative to n_ref
        E = E .* exp( 1i * k0 * (n_real - n_ref) * dz );

        % absorption from extinction coefficient kext
        E = E .* exp( -k0 * n_imag * dz );

        % ---- half diffraction (CN) ----
        E = Aimp \ (Aexp * E);

        E_rz(:,iz) = E;
    end

    I_rz = 0.5 * sp.n * sp.eps0 * sp.c0 * abs(E_rz).^2;

end
