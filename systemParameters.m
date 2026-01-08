function sp = systemParameters()

    sp.c0 = 299792458;  % Speed of light in vacuum [m/s]
    sp.h = 6.62607e-34; % Planck's constant [J*s]
    sp.T = 300;         % Room temperature [K]

    sp.NA = 0.4;
    sp.n = 3.5;         % Refractive index in Silicon
    
    % Sample dimensions:
    sp.Lx = 1e-5;       % 10 micron, [m]
    sp.Lz = 1e-4;       % 100 micron, [m]

    % Diffusion parameters:
    sp.D = 2.6e-4;      % Diffusion parameter [m^2/s]
    sp.tau = 2.5e-9;    % Effective recombination time of 25[ns], [s]


    % Laser parameters:
    sp.wl1 = 1550e-9;   % wavelength [m]
    sp.wl2 = 775e-9;    % wavelength [m]
    sp.D1 = 1.2e-3;     % beam width [?] ??????????
    sp.D2 = 3e-5;       % beam width [?] ??????????
    sp.pulse_energies = (0:0.25:5) .* 1e-6; % Pulse energy (intensity) [J]
    sp.alpha775 = 1e5;  % α for 775nm in n-type Si acc. to BS [1/m]
    sp.alpha1550 = 1e5; % α for 1550nm in n-type Si acc. to BS [1/m] ????????????

    % Penetration depth:
    sp.penDepth = 1e-5; % 10 micron, [m]

end