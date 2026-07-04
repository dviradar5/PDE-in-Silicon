%% CHECKKKK

% =========================================================================
% Fiber Coupling Simulation (1D Continuous GRIN Waveguide)
% Method: Paraxial Split-Step Fourier Method (SSFM)
% =========================================================================
clear; clc; %close all;

%% 1. Tunable Parameters
wl = 0.755e-6;                      % Probe Wavelength [m]

% GRIN & Fiber Profile Parameters
pump_w0 = 4.24e-6;                 % Pump beam waist [m]
l = 1;                              % OAM mode index (LG01)

% Effective Core defined by pump peaks: 2 * sqrt(|l|/2) * w0
effective_core_width = 2 * sqrt(abs(l)/2) * pump_w0; 

% Continuous Background Properties (Bulk Medium)
n_base = 3.7139;                    % Unperturbed background refractive index
delta_n = -0.59;                    % Max index change at the donut peaks
alpha_base_cm = 1296.2;             % Background absorption [1/cm]
delta_alpha_cm = 600;                % Max absorption increase [1/cm]

beam_diameter = 2*2.267e-6;          % Incident probe beam diameter (2*w0) [m]
Lx = 20e-6;                         % Transverse computational domain width [m]
Lz = 25e-6;                         % Total propagation distance [m]
Nz = 500;                           % Number of z-steps

%% 2. Grid Setup
Nx = 2048;                          % High resolution
x = linspace(-Lx/2, Lx/2, Nx);      % Transverse axis [m]
dx = x(2) - x(1);
z = linspace(0, Lz, Nz);            % Propagation axis [m]
dz = z(2) - z(1);

k0 = 2*pi / wl;                     % Vacuum wavenumber

%% 3. Complex Refractive Index Profile (Continuous GRIN)
% LG01 Transverse Shape Function (Normalized to peak at 1)
u = 2 * (x.^2) / (pump_w0^2);
S_x = exp(1) .* u .* exp(-u); 

% Continuous arrays spanning the entire domain
n_array = n_base + (delta_n .* S_x);
alpha_array_cm = alpha_base_cm + (delta_alpha_cm .* S_x);
alpha_array = alpha_array_cm * 100; % Convert [1/cm] to [1/m]

n_ref = n_base;                     % Paraxial reference index 

%% 4. Initial Field (Gaussian Probe)
w0 = beam_diameter / 2;
E = exp(-(x.^2) / w0^2);            % Incident electric field amplitude

%% 5. Split-Step Fourier Pre-allocations
kx = (2*pi/Lx) * [0:(Nx/2-1), (-Nx/2):-1]; 

% Phase operators
diffraction_phase = exp(-1i * (kx.^2) * dz / (2 * k0 * n_ref));
refraction_phase = exp(1i * k0 * (n_array - n_ref) * dz);

% Absorbing Boundary Condition (to stop FFT wrap-around at grid edges)
mask = exp(- (x / (Lx*0.45)).^20 ); 

%% 6. Propagation Loop
I = zeros(Nz, Nx);

for step = 1:Nz
    % Spatial domain (Refraction & localized absorption)
    E = E .* refraction_phase .* mask .* exp(-alpha_array .* dz / 2);
    
    % Fourier domain (Diffraction)
    E_kx = fft(E);
    E_kx = E_kx .* diffraction_phase;
    E = ifft(E_kx);
    
    I(step, :) = abs(E).^2;
end

%% 7. Visualization (Updated with Shaded Core)
% ---------------------------------------------------------
% Figure 1: Verify the Continuous GRIN Profile with Shaded Core
% ---------------------------------------------------------
figure('Color', 'w', 'Position', [50, 50, 700, 400]);

% 1. Plot the data first to establish axis limits
yyaxis left
plot(x*1e6, n_array, 'm-', 'LineWidth', 2.5);
ylabel('Refractive Index', 'FontSize', 12, 'Color', 'm');
set(gca, 'YColor', 'm');

yyaxis right
plot(x*1e6, alpha_array_cm, 'r--', 'LineWidth', 2.5);
ylabel('Absorption [1/cm]', 'FontSize', 12, 'Color', 'r');
set(gca, 'YColor', 'r');

% 2. Retrieve the current axis limits
yLims = ylim; 

% 3. Calculate core bounds
core_limit = effective_core_width / 2;

% 4. Draw the shading using the full vertical range (yLims)
hold on;
fill([-core_limit, core_limit, core_limit, -core_limit]*1e6, ...
     [yLims(1), yLims(1), yLims(2), yLims(2)], ...
     'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Fiber Core');

xlabel('x [\mum]', 'FontSize', 12);
title('Pump-Induced GRIN & Absorption Profile', 'FontSize', 14);
grid on; xlim([-1e6*Lx/2 1e6*Lx/2]);
legend('Location', 'northeast');

% ---------------------------------------------------------
% Figure 2: Beam Evolution
% ---------------------------------------------------------
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

subplot(2,1,1);
imagesc(z*1e6, x*1e6, I');
colormap('hot'); colorbar;
shading interp;
ylabel('x [\mum]', 'FontSize', 12);
xlabel('z [\mum]', 'FontSize', 12);
title('Probe Propagation in GRIN Waveguide', 'FontSize', 14);
hold on;
% Plotting the effective core width lines for reference
plot([0 Lz*1e6], [effective_core_width/2*1e6 effective_core_width/2*1e6], 'w--', 'LineWidth', 1.5);
plot([0 Lz*1e6], [-effective_core_width/2*1e6 -effective_core_width/2*1e6], 'w--', 'LineWidth', 1.5);
ylim([-Lx/2*1e6, Lx/2*1e6]);

subplot(2,1,2);
plot(x*1e6, I(1,:) / max(I(1,:)), 'b-', 'LineWidth', 2);
hold on;
plot(x*1e6, I(end,:) / max(I(end,:)), 'r-', 'LineWidth', 2);
patch([-effective_core_width/2 -effective_core_width/2 effective_core_width/2 effective_core_width/2]*1e6, ...
      [0 1.1 1.1 0], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
legend('Input Probe Intensity', 'Output Probe Intensity', 'Effective Core', 'Location', 'best');
xlabel('x [\mum]', 'FontSize', 12);
ylabel('Normalized Intensity [au]', 'FontSize', 12);
title('Transverse Profile: Entrance vs. Output', 'FontSize', 14);
grid on; axis tight; ylim([0 1.1]); xlim([-Lx/2*1e6, Lx/2*1e6]);