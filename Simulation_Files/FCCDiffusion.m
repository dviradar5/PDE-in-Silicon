%% Explain generation term

function p = FCCDiffusion(pump, t, r, z)
    % FCC diffusion calculator
    % ---------------------------------------------------------------------
    % Numerically calculates spatial and temporal free-charge carriers (FCC)
    % distribution inside the Silicon sample using the diffusion equation
    % in cylindrical coordinates
    %
    % The diffusion equation we solve:
    %                  ∂p/∂t = D*∇²p - p/τ + G(r,z,t)
    % the first term is diffusion, the second is recombination and the
    % third is generation, taking into account FCC's generation time (about
    % 100[ps] according to empirical results) so it won't be instantaneous
    %
    % Uses Crank-Nicolson method in time with operator splitting and the
    % appropriate laplacians with Neumann BC (du/dx=0)
    % =====================================================================
    % INPUTS:
    %        pump - pump laser beam, Laser-type object
    %        t - time vector [s]
    %        r - radial spatial coordinate vector [m] 
    %        z coordinate, propagation vector [m]
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

    % FCC distribution matrix:
    p = zeros(Nr,Nz,Nt);
    
    % The total FCC density distribution that would be created by the pump pulse:
    I = pump.intensityProfileBLDumped(z);     % [W/m^2]
    p_gen = intensityToCarriers(I, pump.pulse_width, pump.lambda); % [m^-3]
    
    % Generation term:
    tau_gen = 100e-12;          % 100[ps] FCC generation time
    t0 = 0;
    g = exp(-4*log(2) * ((t - t0)/tau_gen).^2);
    g_int = trapz(t, g);    % Integrated G
    % w distributes the density to time rations based on the gaussian pulse:
    w = g ./ g_int;         % Normalization

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
        temp1 = Ar \ (p(:,:,i-1) * Bz);     % Nr x Nz
        temp2 = (Br * temp1) / Az;          % Nr x Nz
        
        % Generation:
        temp2 = temp2 + p_gen * (w(i) * dt);
  
        % Recombination:
        %p(:,:,i) = temp2 * (1 - dt/sp.tau);
        p(:,:,i) = temp2 * exp(-dt/sp.tau);
    end

end