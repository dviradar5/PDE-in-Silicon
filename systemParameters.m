%% FINISHED

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
        e = 1.60217663e-19;         % Elementary charge [C]
        me0 = 9.1093837015e-31;     % Electron mass [kg]
        eps0 = 8.854187e-12;        % Vacuum permittivity [F/m]
        mu0 = 1.256637e-6;          % Vacuum permeability [N/A^2]
        T = 300;                    % Room temperature [K]
    
        % Silicon parameters:
        n = 3.7139;                 % Silicon's refractive index at 775[nm]
        ni = 1e10;                  % Intrinsic Si carrier concentration [1/cm^3]
        
        % Penetration depth:
        penDepth = 1e-5;            % 10 micron, [m]
    
        % Sample dimensions:
        Lx = 10e-6;                 % 10 micron [m]
        Lz1 = 20e-6;                % 20 micron [m]
        Lz2 = 25e-6;                % 25 micron [m]
    
        % Diffusion parameters:
        D = 2.6e-4;                 % Diffusion parameter [m^2/s]
        tau = 2.5e-9;               % Effective recombination time [s]
        elt = 50e-12;               % Hot electron life time [s]
    
        % Drude model parameters:
        tau_D = 1e-14;              % Drude damping time [s]
        me_eff = 2.36843976e-31;    % Effective electron mass, 0.26me0 [kg]
        mh_eff = 3.55265964359e-31; % Effective hole mass, 0.39me0 [kg]
        mu_e = 0.14;                % Electron mobility [m^2/(V*s)]
        mu_h = 0.045;               % Hole mobility [m^2/(V*s)]
    
        % Laser and optics parameters:
        NA = 0.4;                   % Objective lens maximal NA
        Dobj = 8e-2;                % Objective lens diameter, [m]
        f = 10e-3;                  % Objective lens focal length, [m]
        wl1 = 1550e-9;              % Wavelength [m]
        wl2 = 775e-9;               % Wavelength [m]
        alpha = 129620;             % α in Silicon at 775[nm] [1/m]
        alpha1550 = 1e5;            % α for 1550nm in n-type Si acc. to BS [1/m]
    
        % Pretty colors for plots:
        colors = ["#000000", "#000080", "#006400", "#1F77B4", "#FF6187",...
                  "#00bfff", "#00A896", "#80B3FF", "#483d8b", "#7B2CBF",...
                  "#9B5DE5", "#8b008b", "#C44569", "#80B3FF", "#696969",...
                  "#8b0000", "#ff4500", "#CC4F1B", "#2f4f4f", "#663399",...
                  "#b8860b", "#E9A820", "#ffd700", "#808000", "#9acd32",...
                  "#c71585", "#E76F51"];
    end

    methods
        function color = randomColor(this)
            % Returns one random color from the color palette
            color = this.colors(randi(numel(this.colors)));
        end
    end
end