function [NAeff, w0] = effectiveNA_from_BE(lambda, NAobj, D0, BE, Dpupil)
    % lambda  - wavelength in vacuum [m]
    % NAobj   - nominal objective NA
    % D0      - original beam diameter before BE [m]
    % BE      - beam-expander magnification, e.g. 1, 1/3, 1/5, 1/10
    % Dpupil  - objective back-pupil diameter [m]

    Dbeam = D0 * BE;

    fillFactor = min(1, Dbeam / Dpupil);

    NAeff = NAobj * fillFactor;

    % Gaussian 1/e^2 waist approximation
    w0 = lambda / (pi * NAeff);
end