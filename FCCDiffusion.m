%% FINISHED

function p = FCCDiffusion(pump, t, r, z)
    % PDE calculation function
    % ---------------------------------------------------------------------
    % Numerically calculates spacial and temporal free-charge carriers (FCC)
    % distribution inside the Silicon sample using diffusion equation
    % =====================================================================
    % INPUT:
    %        pump - pump laser beam. Laser-type object
    %        t - time vector
    %        r - radial spacial coordinate vector [m] 
    %        z - vertical spacial coordinate vector coordiante [m]
    % 
    % OUTPUT:
    %        p - FCC distribution, Nr x Nz x Nt, [1/m^3]
    % *********************************************************************
    sp = systemParameters();

    Nr = length(r);                      % Number of elements
    Nz = length(z);
    Nt = length(t);
    
    % Steps:
    dr = r(2) - r(1);
    dz = z(2) - z(1);
    dt = t(2) - t(1);

    % Charge distribution matrix:
    p = zeros(Nr,Nz,Nt);
    
    % The initial distribution p(r,z,0) relates exactly to the intensity
    % profile of the pump since each photon creates exactly 1 e-h pair:
    I = pump.intensityProfileBLDumped(z);
    p(:,:,1) = intensityToCarriers(I, sp.alpha, pump.lambda) * 1e6; % [m^-3]

    % Crank-Nicolson method with operator splitting:
    Sr = lap1dNeumannCylR(r, dr);           % Nr x Nr
    Sz = lap1dNeumann(Nz, dz);              % Nz x Nz
    
    Ir = speye(Nr);
    Iz = speye(Nz);
    
    % Total laplacian acting on p(r,t):
    Ar = (Ir - 0.5*sp.D*dt*Sr);
    Br = (Ir + 0.5*sp.D*dt*Sr);

    Az = (Iz - 0.5*sp.D*dt*Sz);
    Bz = (Iz + 0.5*sp.D*dt*Sz);

    
    for i = 2:Nt
        % Diffusion:
        temp2 = Ar \ (p(:,:,i-1) * Bz);     % Nr x Nz
        temp3 = (Br * temp2) / Az;          % Nr x Nz
            
        % Recombination:
        p(:,:,i) = temp3 * (1 - dt/sp.tau);
    end
end