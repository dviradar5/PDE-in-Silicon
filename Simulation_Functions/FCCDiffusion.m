%% FINISHED

function p = FCCDiffusion(pump, t, r, z)
    % FCC diffusion calculator
    % ---------------------------------------------------------------------
    % Numerically calculates spatial and temporal free-charge carriers
    % (FCC) distribution inside the Silicon sample using the diffusion
    % equation in cylindrical coordinates
    %
    % The diffusion equation we solve:
    %                  ∂p/∂t = D*∇²p - p/τ + G(r,z,t)
    % the first term is diffusion, the second is recombination and the
    % third is generation, taking into account FCC's generation time (about
    % 100[ps] according to empirical results) so it won't be instantaneous
    %
    % Uses split Crank-Nicolson method in time with operator splitting and
    % the appropriate laplacians with Neumann BC (du/dx=0)
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
    I = pump.intensityProfileBLDumped(z);       % [W/m^2]
    p_gen = intensityToCarriers(I, pump.pulse_width, pump.lambda); % [m^-3]
    
    % Generation term:
    tau_gen = 100e-12;%pump.pulse_width;        % 100[ps] FCC generation time
    t0 = 0;
    
    % Taking the total carriers density created by the entire pulse and
    % distrbuting it in time according to a Gaussian envelope:
    g = exp(-4*log(2) * ((t - t0)/tau_gen).^2);
    g_int = trapz(t, g);    % Integrated generation term G
    w = g ./ g_int;         % Normalization

    % For instant generation:
    %p(:,:,1) = p_gen;

    % Split Crank-Nicolson method with operator splitting:
    Lr = lap1dNeumannCylR(r, dr);                       % Nr x Nr
    Lz = lap1dNeumann(Nz, dz);                          % Nz x Nz
    
    Ir = speye(Nr);
    Iz = speye(Nz);
    
    % Total laplacian acting on p(r,t):
    Mplus_r = (Ir - 0.5*sp.D*dt*Lr);
    Mminus_r = (Ir + 0.5*sp.D*dt*Lr);

    Mplus_z = (Iz - 0.5*sp.D*dt*Lz);
    Mminus_z = (Iz + 0.5*sp.D*dt*Lz);

    for i = 2:Nt
        % Diffusion:
        temp1 = Mplus_r \ (p(:,:,i-1) * Mminus_z);      % Nr x Nz
        temp2 = (Mminus_r * temp1) / Mplus_z;           % Nr x Nz
        
        % Generation (not instant):
        temp2 = temp2 + p_gen * (w(i) * dt);
  
        % Recombination:
        %p(:,:,i) = temp2 * (1 - dt/sp.tau);
        p(:,:,i) = temp2 * exp(-dt/sp.tau);             % Exact solution
    end
end

%% Helper:
function N = intensityToCarriers(I, tau, lambda)
    % Converts intensity to carrier concentration
    % ---------------------------------------------------------------------
    % Calculates electron/hole concentration given the intensity profile
    % 
    % Assumes quantum efficiency of 1 (also means that Ne = Nh)
    % 
    % Fluence calculation takes into account the entire pulse:
    %              I(r,z,t) = I(r,z) * exp(−4ln(2)*(t/τ)^2​)
    %         => F = ∫​I(r,z,t)dt = I(r,z) * sqrt(pi/ln(2)) * τ/2
    % =====================================================================
    % INPUTS:
    %        I - intensity spatial profile, Nr x Nz
    %        tau - gaussian envelope duration [s] 
    %        lambda - beam's wavelength [m]
    % OUTPUT:
    %        N - carrier concentration [m^-3]
    % *********************************************************************

    sp = systemParameters();
    
    F = I * (tau/2) *sqrt(pi/log(2));       % [J/m^2]

    Eph = sp.h * sp.c0 / lambda;                % [J]

    N = sp.alpha .* F / Eph;                % [1/m^3]

end