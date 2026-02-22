%% CHECK SIZES AND UNITS

classdef systemParameters
    % Parameters structure
    % ---------------------------------------------------------------------
    % Keeps all the constants and parameters needed for the simulation in
    % one place, (most) in SI units
    % *********************************************************************
    
    properties
        % Natural constants and parameters:
        c0 = 299792458;             % Speed of light in vacuum [m/s]
        h = 6.62607e-34;            % Planck's constant [J*s]
        kB = 1.380649e-23;          % Boltzmann's constant [J/K]
        e_m0 = 9.1093837015e-31;    % Electron mass [kg]
        eps0 = 8.854187e-12;        % Vacuum permittivity [F/m]
        mu0 = 1.256637e-6;          % Vacuum permeability [N/A^2]
        T = 300;                    % Room temperature [K]
    
        % Silicon parameters:
        n = 3.7139;                 % Silicon's refractive index at 775[nm]
        ni = 1e10;                  % Intrinsic Si carrier concentration [1/cm^3]
        
        % Penetration depth:
        penDepth = 1e-5;            % 10 micron, [m]
    
        % Sample dimensions:
        Lx = 1e-5;                  % 10 micron, [m]
        Lz = 2.5e-5;                % 25 micron, [m]
    
        % Diffusion parameters:
        D = 2.6e-4;                 % Diffusion parameter [m^2/s]
        tau = 2.5e-9;               % Effective recombination time [s]
        elt = 50e-12;               % Hot electron life time [s]
    
        % Drude model parameters:
        tau_D = 1e-14;              % Drude damping time [s]
        mu_e = 1400e-4;             % Electron mobility m^2/(V*s)
        mu_h = 450e-4;              % Hole mobility m^2/(V*s)
    
        % Laser parameters:
        NA = 0.4;
        wl1 = 1550e-9;              % Wavelength [m]
        wl2 = 775e-9;               % Wavelength [m]
        D1 = 1.2e-3;                % Beam width [?] ??????????
        D2 = 3e-5;                  % Beam width [?] ??????????
        alpha = 129620;             % α in Silicon at 775[nm] [1/m]
        %alpha1550 = 1e5;            % α for 1550nm in n-type Si acc. to BS [1/m]
    end
end