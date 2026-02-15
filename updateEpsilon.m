%% CHANGE THIS
%% BASED ON ....
%% DOCUMENT

function [eps_rz, n_complex, dn, dalpha] = updateEpsilon(Ne, Nh, lambda, n_base, alpha_base)
    % updateEpsilon  Update complex permittivity from free-carrier densities.
    %
    % Uses your BennetSoref(Ne,Nh,lambda) to get:
    %   dn      : change in refractive index (dimensionless)
    %   dalpha  : change in absorption coefficient [1/m]
    %
    % Then converts absorption coefficient to extinction coefficient:
    %   k = alpha * lambda / (4*pi)
    %
    % And returns:
    %   n_complex = (n_base + dn) - 1i*(k_base + dk)
    %   eps_rz    = n_complex.^2
    %
    % Inputs:
    %   Ne, Nh      carrier densities [cm^-3] (Nr x Nz) or any grid
    %   lambda      wavelength [m]
    %   n_base      base refractive index (e.g. 3.48 at 1550nm)
    %   alpha_base  base absorption coefficient [1/m] (can be 0 if negligible)
    %
    % Outputs:
    %   eps_rz      complex permittivity
    %   n_complex   complex refractive index
    %   dn          delta n
    %   dalpha      delta alpha [1/m]
    % *********************************************************************

    arguments
        Ne
        Nh
        lambda (1,1) double {mustBePositive}
        n_base (1,1) double {mustBePositive}
        alpha_base (1,1) double {mustBeNonnegative}
    end

    % Your empirical carrier-induced effects:
    [dn, dalpha] = BennetSoref(Ne, Nh, lambda);  % dn dimensionless, dalpha [1/m]

    % Convert absorption coefficient -> extinction coefficient k
    k_base = alpha_base * lambda / (4*pi);
    dk     = dalpha     * lambda / (4*pi);

    % Complex refractive index and epsilon:
    n_complex = (n_base + dn) - 1i*(k_base + dk);
    eps_rz    = n_complex.^2;   % mu_r = 1
end
