%% CHECK BL DUMPING

classdef laser
    % Laser Beam Class
    % ---------------------------------------------------------------------
    % Describes single laser beam with spatial and temporal Gaussian pulse
    % *********************************************************************

    properties
        lambda;              % Wavelength [m]
        f;                   % Frequency [Hz]
        ang_f;               % Angular frequency [rad/s]

        pulse_width;         % Pulse width (intensity FWHM) [s]
        pulse_energy;        % Pulse total energy [J]

        %rep_rate;            % Repetition rate [Hz]
        %t0;                  % Pulse center time [s]

        %beam_D;              % Beam diameter [m]

        type;                % Gaussian or "Donut"
        w0;                  % Waist radius at z=z0 [m]
        z0;                  % Waist location along z axis [m]
        E0;                  % Peak field amplitude at waist center [V/m]

        profile;             % Spatial field profile matrix E(r,z) [V/m]
    end

    methods
        % Constructor:
        function this = laser(wl, pls_wd, pls_energy, type, r, phi, z, w0, z0)
            % Constructs a laser beam
            % =============================================================
            % INPUTS:
            %        wl - wavelength [m]
            %        pls_wd - pulse width (intensity FWHM) [s]
            %        pls_energy - pulse total energy [J]
            %        type - beam spatial profile type, "Gauss" or "Donut"
            %        r - radial coordinate vector [m]
            %        phi - azimutal coordinate vector [rad]
            %        z -z coordinate, propagation vector [m]
            %        w0 - waist radius at z = z0 [m]
            %        z0 - waist location along z axis [m]
            % *************************************************************

            sp = systemParameters();

            this.lambda = wl;
            this.f = sp.c0 / wl;
            this.ang_f = 2*pi*this.f;

            this.pulse_width  = pls_wd;
            this.pulse_energy = pls_energy;

            this.type = string(type);
            this.w0 = w0;
            this.z0 = z0;

            % Computing correct peak field amplitude E0:
            this.E0 = this.computeE0FromEnergy();

            % Building spatial field profile:
            this.profile = this.beamProfile(this.type, r, phi, z, this.w0, this.z0, this.E0);
        end
        
        function E0 = computeE0FromEnergy(this)
            % Computes peak electric field amplitude
            % -------------------------------------------------------------
            % Calculates E0, the electric field's amplitude at the pulse's
            % peak at the waist given the pulse's total energy
            % 
            % The intensity spatial profile at the waist:
            %               I(r) = I0 * exp(-2*(r/w0)^2)
            %              I0 = 0.5 * n * eps0 * c * |E0|^2
            %
            % The intensity temporal Gaussian profile is:
            %                g(t) = exp(-4ln(2)*(t/τ)^2)
            % where τ is intensity FWHM
            %
            % The total pulse's energy is given by integrating I(r,t) in
            % both space and time:
            %           E = I0 * (pi*τ*w0^2/4) * sqrt(pi/ln(2))
            % =============================================================
            % INPUT:
            %        this - this laser-type object
            % OUTPUT:
            %        E0 - peak field amplitude at waist center [V/m]
            % *************************************************************
            
            sp = systemParameters();
            
            if this.pulse_energy == 0
                E0 = 0;
                return;
            end
            
            % Results of the integrals:
            spatialInt = pi*(this.w0^2)/2;
            temporalInt = this.pulse_width*sqrt(pi/log(2))/2;

            I0 = this.pulse_energy / (spatialInt * temporalInt);    % [W/m^2]
            E0 = sqrt(2*I0/(sp.n*sp.eps0*sp.c0));                   % [V/m]
        end
        
        function energy = computeEnergyFromE0(this)
            % Computes pulse's energy
            % -------------------------------------------------------------
            % Calculates the pulse's total energy given E0, the electric
            % field's amplitude peak at the waist
            % 
            % The total pulse's energy is given by integrating I(r,t) in
            % both space and time:
            %           E = I0 * (pi*τ*w0^2/4) * sqrt(pi/ln(2))
            % where
            %              I0 = 0.5 * n * eps0 * c * |E0|^2
            % =============================================================
            % INPUT:
            %        this - this laser-type object
            % OUTPUT:
            %        energy - pulse total energy [J]
            % *************************************************************
            
            sp = systemParameters();
            
            if this.E0 == 0
                energy = 0;
                return;
            end
            
            I0 = 0.5 * sp.n * sp.eps0 * sp.c0 * abs(this.E0)^2; % [W/m^2]
            
            energy = I0 * (pi*this.pulse_width*this.w0^2/4) * sqrt(pi/ln(2));    % [J]
        end

        function prf = beamProfile(this, type, r, phi, z, w0, z0, E0)
            % Beam profile creator
            % -------------------------------------------------------------
            % Creates and returns appropriate spatial field profile E(r,z)
            % =============================================================
            % INPUTS:
            %        this - this laser-type object
            %        type - beam spatial profile type, "Gauss" or "Donut"
            %        r - radial coordinate vector [m]
            %        phi - azimutal ????????????????
            %        z -z coordinate, propagation vector [m]
            %        w0 - waist radius at z=z0 [m]
            %        z0 - waist location along z axis [m]
            %        E0 - peak field amplitude at waist center [V/m]
            % OUTPUT:
            %        prf - spatial profile matrix, Nr x Nz, [V/m]
            % *************************************************************

            if string(type) == "Gauss"
                prf = GB(r, z, this.lambda, w0, z0, E0);
            
            elseif string(type) == "Donut"
                prf = LGB01(r, phi, z, this.lambda, w0, z0, E0);
            
            else
                error("Inappropriate type. Please use 'Gauss' or 'Donut'");
            end
        end
        
        function I = intensityProfileBLDumped(this, z)  % Maybe change the alpha when it changes
            % Builds intensity matrix including Beer-Lmbert's law
            % -------------------------------------------------------------
            % Creates and returns peak spatial intensity including
            % Beer-Lambert's law:
            %             Imaterial(r,z) = I(r,z)*exp(-α*z)
            % =============================================================
            % INPUTS:
            %        this - this laser-type object
            %        z -z coordinate, propagation vector [m]
            % OUTPUT:
            %        I - spatial intensity profile matrix, Nr x Nz, [W/m^2]
            % *************************************************************
            
            sp = systemParameters();
            
            z = z(:).';     % 1 x Nz
    
            I = 0.5 * sp.n * sp.eps0 * sp.c0 * abs(this.profile).^2 .* exp(-sp.alpha * z);
            %I = I / max(I(:));
        end
        
        % Update function:
        function this = updateFieldProfile(this, r, phi, z, w0, z0, E0, type, pls_wd, pulse_energy)
            % Updates electric field spatial profile (and other parameters)
            % =============================================================
            % INPUTS:
            %        this - this laser-type object
            %        r - radial coordinate vector [m]
            %        phi - azimutal coordinate vector [rad]
            %        z -z coordinate, propagation vector [m]
            %        w0 - waist radius at z = z0 [m]
            %        z0 - waist location along z axis [m]
            %        E0 - peak field amplitude at waist center [V/m]
            %        type - beam spatial profile type, "Gauss" or "Donut"
            %        pls_wd - pulse width (intensity FWHM) [s]
            %        pulse_energy - pulse total energy [J]
            % OUTPUT:
            %        this - this updated laser-type object
            % *************************************************************

            arguments
                this
                r
                phi
                z
                w0 double = []              % optional
                z0 double = []              % optional
                E0 double = []              % optional (if provided, overrides energy/tau)
                type string = ""            % optional
                pls_wd double = []             % optional pulse width (intensity FWHM) [s]
                pulse_energy double = []    % optional pulse energy [J]
            end

            w0_changed = ~isempty(w0);
            z0_changed = ~isempty(z0);
            E0_changed = ~isempty(E0);
            type_changed = strlength(type) > 0;
            width_changed = ~isempty(pls_wd);
            Energy_changed = ~isempty(pulse_energy);
        
            % Update stored parameters when provided:
            if z0_changed, this.z0 = z0;
            else, z0 = this.z0; end
        
            if w0_changed, this.w0 = w0;
            else, w0 = this.w0; end
        
            if width_changed, this.pulse_width = pls_wd;
            else, pls_wd = this.pulse_width; end
        
            if Energy_changed, this.pulse_energy = pulse_energy;
            else, pulse_energy = this.pulse_energy; end
        
            if type_changed, this.type = type;
            else, type = this.type; end
        
            % E0 depends on other parameters:
            if E0_changed
                this.E0 = E0;
        
                % Update pulse energy based on the new E0, w0 and width:
                this.pulse_energy = this.computeEnergyFromE0();
        
            else
                % If energy/width/w0 changed, E0 changes too:
                if Energy_changed || width_changed || w0_changed
                    this.E0 = this.computeE0FromEnergy();
                end
                
                E0 = this.E0;
            end
        
            % Updating the profile:
            this.profile = beamProfile(type, r, phi, z, w0, z0, E0);
        end
    end
end
