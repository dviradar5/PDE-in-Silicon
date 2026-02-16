%% CHECK SIZES AND UNITS

function sp = systemParameters()
    % Keeps all the constants and parameters needed for the simulation in 
    % one place
    % *********************************************************************

    % Natural constants and parameters:
    sp.c0 = 299792458;  % Speed of light in vacuum [m/s]
    sp.h = 6.62607e-34; % Planck's constant [J*s]
    sp.kB = 1.380649e-23;       % Boltzmann's constant [J/K]
    sp.e_m0 = 9.1093837015e-31; % Electron mass [kg]
    sp.eps0 = 8.854187e-12; % Vacuum permittivity [F/m]
    sp.mu0 = 1.256637e-6;   % Vacuum permeability [N/A^2]
    sp.T = 300;         % Room temperature [K]

    % Silicon parameters:
    sp.n = 3.48;      % Refractive index in Silicon (3.9766)
    sp.ni = 1e10;       % Intrinsic Si carrier concentration [1/cm^3]
    
    % Penetration depth:
    sp.penDepth = 1e-5; % 10 micron, [m]

    % Sample dimensions:
    sp.Lx = 1e-5;       % 10 micron, [m]
    sp.Lz = 2.5e-5;     % 25 micron, [m]

    % Diffusion parameters:
    sp.D = 2.6e-4;      % Diffusion parameter [m^2/s]
    sp.tau = 2.5e-9;    % Effective recombination time of 25[ns], [s]
    sp.elt = 50e-12;    % hot election life time [s]

    % Drude model parameters:
    sp.tau_D = 1e-14;   % Drude damping time [s]
    sp.mu_e = 1400e-4;  % electron mobility m^2/(V?s) ?????????
    sp.mu_h = 450e-4;   % hole mobility m^2/(V?s) ????????

    % Laser parameters:
    sp.NA = 0.4;
    sp.wl1 = 1550e-9;   % wavelength [m]
    sp.wl2 = 775e-9;    % wavelength [m]
    sp.D1 = 1.2e-3;     % beam width [?] ??????????
    sp.D2 = 3e-5;       % beam width [?] ??????????
    sp.pulse_energies = (0:0.25:5) .* 1e-6; % Pulse energy (intensity) [J]
    sp.alpha = 129620;  % α in Si [1/m] ???????????????????????
    sp.alpha775 = 1e5;  % α for 775nm in n-type Si acc. to BS [1/m] ??????????
    %sp.alpha1550 = 1e5; % α for 1550nm in n-type Si acc. to BS [1/m] ????????????


end