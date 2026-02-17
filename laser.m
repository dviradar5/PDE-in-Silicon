%% ADD PULSE AND TIME DEPENDANCE, ADD ENERGY
%% DOCUMENT
%% CHECK FUNCTIONS

classdef laser
    % Laser Beam Class
    % ---------------------------------------------------------------------
    % Describes single laser beam (pump/probe) with spatial + temporal pulse
    % scaling from pulse energy.
    %
    % Conventions used here:
    %   - GB/LGB01 return complex electric field E(r,z) [V/m] when given E0 [V/m]
    %   - pulseEnergy is energy per pulse [J]
    %   - pulse_width is INTENSITY FWHM duration [s] (common in optics)
    %   - temporal envelope is Gaussian in time
    %
    % *********************************************************************

    properties
        lambda;              % Wavelength [m]
        f;                   % Frequency [Hz]
        ang_f;               % Angular frequency [rad/s]

        pulse_width;         % Intensity FWHM [s]
        pulse_energy;        % Pulse energy [J]

        rep_rate;            % Repetition rate [Hz] (optional; default 0 => single pulse)
        t0;                  % Pulse center time [s] (default 0)

        beam_D;              % Beam diameter [m] (optional, not used directly)

        type;                % "Gauss" or "Donut"
        w0;                  % Waist radius at z=z0 [m]
        z0;                  % Waist location along z [m]


        E0;                  % Peak field amplitude at waist center [V/m]

        profile;             % Spatial complex field profile E(r,z) [V/m] (for one pulse peak)
    end

    methods
        function obj = laser(wl, pls_wd, pls_enr, bm_D, type, r, phi, z, w0, z0, varargin)
            % Constructs a laser beam
            %
            % Inputs:
            %   wl        wavelength [m]
            %   pls_wd    pulse width INTENSITY FWHM [s]
            %   pls_enr   pulse energy [J]
            %   bm_D      beam diameter [m] (optional)
            %   type      "Gauss" or "Donut"
            %   r,phi,z   grids
            %   w0,z0     waist and waist location
            %
            % Optional name-value pairs:
            %   'rep_rate'  [Hz] default 0 (single pulse)
            %   't0'        [s]  default 0

            sp = systemParameters();

            obj.lambda = wl;
            obj.f = sp.c0 / wl;
            obj.ang_f = 2*pi*obj.f;

            obj.pulse_width  = pls_wd;
            obj.pulse_energy = pls_enr;

            obj.beam_D = bm_D;

            obj.type = string(type);
            obj.w0 = w0;
            obj.z0 = z0;

            % defaults
            obj.rep_rate = 0;
            obj.t0 = 0;

            % parse optional name-value args
            if ~isempty(varargin)
                for k = 1:2:numel(varargin)
                    name = lower(string(varargin{k}));
                    val  = varargin{k+1};
                    switch name
                        case "rep_rate"
                            obj.rep_rate = val;
                        case "t0"
                            obj.t0 = val;
                        otherwise
                            error("Unknown option '%s'. Valid: 'rep_rate','t0'.", name);
                    end
                end
            end

            % Compute correct E0 from pulse energy, width, waist
            obj.E0 = obj.computeE0FromEnergy();

            % Build spatial field profile at pulse peak (temporal envelope = 1)
            obj.profile = obj.beamProfile(obj.type, r, phi, z, obj.w0, obj.z0, obj.E0);
        end

        function prf = beamProfile(this, type, r, phi, z, w0, z0, E0)
            % Creates and returns appropriate spatial field profile E(r,z) [V/m]
            type = string(type);

            if type == "Gauss"
                prf = GB(r, z, this.lambda, w0, z0, E0);
            elseif type == "Donut"
                prf = LGB01(r, phi, z, this.lambda, w0, z0, E0);
            else
                error("Inappropriate type. Use 'Gauss' or 'Donut'.");
            end
        end

        function E0 = computeE0FromEnergy(this)
            % Compute E0 [V/m] so that the pulse energy matches pulse_energy.
            %
            % Assumptions:
            %   - intensity spatial profile at waist: I(r) = I0 * exp(-2 r^2 / w0^2)
            %   - intensity temporal profile: g_I(t) = exp(-4 ln2 * t^2 / tau^2)
            %   - pulse_width = tau is INTENSITY FWHM
            %
            % Then:
            %   E_pulse = I0 * (pi w0^2 / 2) * (tau/2)*sqrt(pi/ln2)
            %   I0 = 0.5 * n * eps0 * c * |E0|^2
            %
            % IMPORTANT:
            %   This scaling is defined at the WAIST plane.
            %   If z0 is not inside the sample, you're still fine as long as
            %   you treat E0 as waist-center amplitude.
            
            sp = systemParameters();

            tau = this.pulse_width;
            w0  = this.w0;

            if tau <= 0 || w0 <= 0
                error("pulse_width and w0 must be positive.");
            end
            if this.pulse_energy < 0
                error("pulse_energy must be nonnegative.");
            end

            spatialInt  = (pi*w0^2)/2;
            temporalInt = (tau/2)*sqrt(pi/log(2));

            if this.pulse_energy == 0
                E0 = 0;
                return;
            end

            I0 = this.pulse_energy / (spatialInt * temporalInt);  % [W/m^2]
            E0 = sqrt(2*I0/(sp.n*sp.eps0*sp.c0));              % [V/m]
        end

        function a = temporalEnvelopeField(this, t)
            % Field envelope a(t) for a Gaussian pulse or pulse train.
            % Returns dimensionless a(t) such that peak = 1.
            %
            % Field envelope for intensity-FWHM tau:
            %   I(t) ~ exp(-4 ln2 (t/tau)^2)
            % so field envelope is sqrt(I):
            %   a(t) ~ exp(-2 ln2 (t/tau)^2)

            tau = this.pulse_width;
            if tau <= 0
                error("pulse_width must be positive.");
            end

            if this.rep_rate <= 0
                % single pulse
                a = exp(-2*log(2)*((t - this.t0).^2)/(tau^2));
                return;
            end

            % pulse train
            Trep = 1/this.rep_rate;
            a = zeros(size(t));

            tmin = min(t); tmax = max(t);
            m_min = floor((tmin - this.t0)/Trep) - 2;
            m_max = ceil( (tmax - this.t0)/Trep) + 2;

            for m = m_min:m_max
                tc = this.t0 + m*Trep;
                a = a + exp(-2*log(2)*((t - tc).^2)/(tau^2));
            end
        end
        
        function E = updateProfile(this, profile)
            this.profile = profile;
            E = this.profile;
        end

        function E = fieldProfile(this, t)
            % Returns full spatio-temporal field E(r,z,t):
            %   E(r,z,t) = profile(r,z) * a(t)
            %
            % Output: Nr x Nz x Nt
            a = this.temporalEnvelopeField(t); % 1 x Nt (or Nt x 1)
            E = this.profile .* reshape(a, 1, 1, []);
        end

        function I = intensityProfilePhysical(this)
            % Physical intensity at pulse peak (a(t)=1):
            % I = 0.5 n eps0 c |E|^2  [W/m^2]
            eps0 = 8.854187817e-12;
            c0   = 299792458;
            I = 0.5*this.n_medium*eps0*c0 * abs(this.profile).^2;
        end

        function I = intensityProfilePhysical_t(this, t)
            % Physical intensity I(r,z,t) [W/m^2], Nr x Nz x Nt
            eps0 = 8.854187817e-12;
            c0   = 299792458;

            a = this.temporalEnvelopeField(t);
            I_rz = 0.5*this.n_medium*eps0*c0 * abs(this.profile).^2; % Nr x Nz
            I = I_rz .* reshape(abs(a).^2, 1, 1, []);
        end

        function [I, Imax] = intensityProfileNormalized(this)
            % Normalized spatial intensity (peak = 1), Nr x Nz
            I = abs(this.profile).^2;
            Imax = max(I(:));
            I = I / Imax;
        end

        function I = intensityProfileBLDumped(this, z)
            % Spatial intensity including Beer-Lambert attenuation:
            % I(r,z) = I(r,z)*exp(-alpha_base*z)
            %
            % Inputs:
            %   zvec       z vector [m]
            %   alpha_base absorption coefficient [1/m]
            %
            % NOTE: Using exp(-alpha*z) assumes homogeneous alpha_base.
            sp = systemParameters();
            
            z = z(:).'; % 1 x Nz
            I0 = 0.5 * sp.n * sp.eps0 * sp.c0 * abs(this.E0)^2;   % [W/m^2]
    
            % Spatial factor from field profile (already contains E0), but use |E|^2 / |E0|^2
            spatial = abs(this.profile).^2 / (abs(this.E0)^2);     % dimensionless
    
            I = I0 * spatial .* exp(-sp.alpha * z);             % [W/m^2] peak (t at pulse center)
            %I = abs(this.profile).^2 .* exp(-sp.alpha * z);
            %I = I / max(I(:));
        end
        
        % function I = intensityXY(Ir, r, x)
        %     [X,Y] = meshgrid(x,x);
        %     R = hypot(X,Y);
        %     I = interp1(r, Ir, R, 'linear', 0);
        % end

        function Trep = period(this)
            if this.rep_rate <= 0
                Trep = Inf;
            else
                Trep = 1/this.rep_rate;
            end
        end
    end
end
