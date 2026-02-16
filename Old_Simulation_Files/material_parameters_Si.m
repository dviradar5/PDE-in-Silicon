m_ce = 0.26 * 9.11e-31; % Example electron effective mass
m_ch = 0.38 * 9.11e-31; % Example hole effective mass
wvlen = linspace(0.4e-6, 1.0e-6, 100); % Wavelength range
n = 3.5 * ones(size(wvlen)); % Example refractive index
k = 0.02 * ones(size(wvlen)); % Example extinction coefficient
eps = (n - 1i * k).^2; % Compute permittivity

save('material_parameters_Si.mat', 'm_ce', 'm_ch', 'wvlen', 'n', 'k', 'eps');
%change 2025