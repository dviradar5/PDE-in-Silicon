function lambda = getProbeWavelength(probe)
%GETPROBEWAVELENGTH Robustly extract wavelength [m] from a probe struct/object.
%
% Supported field/property names: lambda, wl, wl0, wavelength, wl2.
% The function intentionally errors if no wavelength is found, because using
% a hidden default can silently change diffraction/focusing results.

    candidates = {'lambda','wl','wl0','wavelength','wl2'};
    lambda = [];

    for k = 1:numel(candidates)
        name = candidates{k};
        if isstruct(probe) && isfield(probe, name)
            lambda = probe.(name);
            break;
        end
        if isobject(probe) && isprop(probe, name)
            lambda = probe.(name);
            break;
        end
    end

    if isempty(lambda)
        error(['Could not find probe wavelength. Add one of these fields/properties ', ...
               'to probe: lambda, wl, wl0, wavelength, wl2. Units must be meters.']);
    end

    lambda = double(lambda);
    if ~isscalar(lambda) || ~isfinite(lambda) || lambda <= 0
        error('Probe wavelength must be a positive finite scalar in meters.');
    end

    % Safety check only: most optical wavelengths in this project are < 10 um.
    if lambda > 1e-3
        warning(['Probe wavelength is larger than 1 mm. I assume it is already in meters. ', ...
                 'If you supplied nm or um, convert it before calling the solver.']);
    end
end
