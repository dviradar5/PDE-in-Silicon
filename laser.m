%% ADD PULSE AND TIME DEPENDANCE

classdef laser
    % Laser Beam Class
    % ---------------------------------------------------------------------
    % Describes single laser beam (i.e., pump/probe) to all its properties
    % and allows parameter manipulations and other initializations
    % *********************************************************************

    properties
        lambda;             % Wavelength [m]
        f;                  % Frequency, [Hz]
        ang_f;              % Anguar frequency [rad/s]
        pulse_width;        % Pulse width [s]
        beam_D;             % Beam diameter [m]         DO SOMETHING
        profile;            % Spacial beam profile
    end

    methods
        function obj = laser(wl, pls_wd, bm_D, type, r, phi, z, w0, z0, E0)
            % Constructs a laser beam
            % *************************************************************
            sp = systemParameters();

            % Updates properties:
            obj.lambda = wl;
            obj.f = sp.c0 / wl;
            obj.ang_f = 2 * pi * obj.f;
            obj.pulse_width = pls_wd;
            obj.beam_D = bm_D;
            obj.profile = beamProfile(obj, type, r, phi, z, w0, z0, E0);
        end

        function prf = beamProfile(obj, type, r, phi, z, w0, z0, E0)
            % Creates and returns apropriate spacial profile in cylindrical
            % coordinates
            % *************************************************************

            if strcmp(type, "Gauss")        % Gaussian beam
                prf = GB(r, z, obj.lambda, w0, z0, E0);

            elseif strcmp(type, "Donut")    % Lagguere-Gauss beam of mode 01  
                prf = LGB01(r, phi, z, obj.lambda, w0, z0, E0);

            else
                error("Inapropriate type; Please try again with a valid type");
            end
        end
    end
end