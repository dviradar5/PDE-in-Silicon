%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN SIMULATION FILE
% -------------------------------------------------------------------------
% Defines all the necessary structures, conducts the experiment flow, plots
% graphs and calculates necessary 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paraneters and Vector Definitions
% -------------------------------------------------------------------------

sp = systemParameters();

% Temporal vector intialization:
tf = 3000;                  % 3[ns]
t = 0:10:tf;
t = t .* 1e-12;             % Time vector, 0-3[ns], [s]
Nt = numel(t);

% Spatial vector intialization:
Nr = 201;                               % Number of elements
Nz = 201;

z = linspace(0, sp.Lz, Nz);

r = linspace(0,sp.Lx,Nr);       % 10x10 micron sample
phi = atan(1);                          % y = x, radial symmetry

x = linspace(-sp.Lx, sp.Lx, 2*Nr-1);
Nx = numel(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pump Beam
% -------------------------------------------------------------------------
% NOTE: pump doesn't need to be t-dependant. Time in this simulation can be
% regarded as the time since when the pump hit the sample.
% =========================================================================

% Pump parameters:
pump_width = 3e-11;     % Pulse of 30[ps]
pump_E = 35e-9;         % Beam total energy [J]
w0 = 5e-6;
z0 = 0;                 % Beam waist location

% Pump laser:
pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, w0, z0);  % 775[nm]
Ipump = pump.intensityProfileBLDumped(z);
Ipump_xz = cylToCart(Ipump,r,x);

% Plotting pump intensity:
%PF_plot_xz(Ipump_xz, z, x, "Pump Intensity [W/m^2]");
%PF_plot_xz(Ipump_xz/max(Ipump_xz(:)), z, x, "Pump Normalized");

iz = 1;     % z=0
%PF_beamIntensityPlot_xy(Ipump(:,iz), r, x, z(iz));
%PF_beamIntensityPlot_xy(Ipump(:,1)/max(Ipump(:,1)), r, x, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FCC Creation & Diffusion
% -------------------------------------------------------------------------

% Shining pump beam on the sample, causing e-h generation and diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)

% Finding the time when we get maximal concentration:
[~, linIdx] = max(pDiff(:));
[~, ~, itMax] = ind2sub(size(pDiff), linIdx);   % itMax = 12 = 110[ps]

%pxz = cylToCart(pDiff(:,:,itMax),r,x) * 1e-6;
%PF_plot_xz(pxz,z,x,"FCC Maximal Concentration in [1/cm^3] at t=110[ps]");


%tIdx = [1, 19, 20, 50, 100, 150, 270, Nt];
%PF_colormapAnimation_xz(pDiff*1e-6, r, z ,x, 1:Nt, "FCC Concentration vs. Time", "FCC concentration [1/cm^3]", false, 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refractive Index Calculation
% -------------------------------------------------------------------------

% Calculating the complex refractive index n(r,z,t) (changes due to FCC generation):
n_complex = complexRefractiveIndex(pDiff, pump.lambda);

%PF_complexRefractiveIndex(n_complex(:,:,itMax), r, z, x, pump.lambda, iz, t(itMax));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probe Beam
% -------------------------------------------------------------------------

% Probe laser:
probe_wd = 5e-11;   % Pulse of 50[ps]
probe_E = pump.pulse_energy/100;        % [J] 
probe_w0 = 2*w0;

probe = laser(sp.wl2, probe_wd, probe_E, "Gauss", r, phi, z, probe_w0, 0);   % 775[nm]
Iprobe = probe.intensityProfileBLDumped(z);
Iprobe_xz = cylToCart(Iprobe,r,x);

%PF_plot_xz(Iprobe_xz, z, x, "Probe Intensity [W/m^2]");
%PF_plot_xz(Iprobe_xz/max(Iprobe_xz(:)), z, x, "Probe Normalized");

%PF_FWHM_z(FWHM(Iprobe,r),z,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probe Propogation (BPM)
% -------------------------------------------------------------------------

Ipropagate = zeros(Nr,Nz,Nt);           % Probe intensity after propogation
Ipropagate_xz = zeros(Nx,Nz,Nt);

for i = 1:Nt
    [~, Ipropagate(:,:,i)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,i));
    Ipropagate_xz(:,:,i) = cylToCart(Ipropagate(:,:,i),r,x);
end

%PF_plot_xz(Ipropagate_xz(:,:,itMax), z, x, "Probe Intensity with 110[ps] Delay at Maximum Concentration");
%PF_plot_xz(Ipropagate_xz/max(Ipropagate_xz(:)), z, x, "Probe Normalized");

% Finding the z index where we get maximal intensity:
slice = Ipropagate_xz(:,:,itMax);
[~, linIdx] = max(slice(:));
[~, izMax] = ind2sub(size(slice), linIdx);

%PF_intensity_x(Ipropagate_xz(:,izMax,itMax),x,"Probe Intensity at Maximum Focusing and 110[ps] delay");
%PF_colormapAnimation_xz(Ipropagate,r,z,x,1:Nt,"Probe Propogation", "Intensity [W/m^2]", false, 50)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons and Other Calculations
% -------------------------------------------------------------------------

% Comparing the probe with and without pump:
%PF_probeComparison(Iprobe_xz, Ipropagate_xz(:,:,itMax), x, z, t(itMax), izMax);

% Checking the maximal intensity vs. delay time:
%PF_maxIntensity(Ipropagate_xz, t);   

% FWHM calculation:
fwhm = zeros(Nz,Nt);

for i = 1:Nt
    fwhm(:,i) = FWHM(Ipropagate(:,:,i),r);
end

PF_FWHM_z(fwhm(:,itMax),z,t(itMax));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pump Energy Related Effects
% -------------------------------------------------------------------------

% Different pump energies:
% pump16 = laser(sp.wl2, pump_width, 16e-9, "Donut", r, phi, z, w0, z0);
% Ipump16_xz = cylToCart(pump16.intensityProfileBLDumped(z),r,x);
% pump64 = laser(sp.wl2, pump_width, 64e-9, "Donut", r, phi, z, w0, z0);
% Ipump64_xz = cylToCart(pump64.intensityProfileBLDumped(z),r,x); 
% pump128 = laser(sp.wl2, pump_width, 128e-9, "Donut", r, phi, z, w0, z0);
% Ipump128_xz = cylToCart(pump128.intensityProfileBLDumped(z),r,x);
% 
% n_complex16 = complexRefractiveIndex(FCCDiffusion(pump16, t, r, z), pump.lambda);
% n_complex64 = complexRefractiveIndex(FCCDiffusion(pump64, t, r, z), pump.lambda);
% n_complex128 = complexRefractiveIndex(FCCDiffusion(pump128, t, r, z), pump.lambda);
% 
% Ipropagate16 = zeros(Nr,Nz,Nt); Ipropagate64 = zeros(Nr,Nz,Nt); Ipropagate128 = zeros(Nr,Nz,Nt);
% Ipropagate16_xz = zeros(Nx,Nz,Nt); Ipropagate64_xz = zeros(Nx,Nz,Nt); Ipropagate128_xz = zeros(Nx,Nz,Nt);
% 
% [~, Ipropagate16(:,:,itMax)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex16(:,:,itMax));
% Ipropagate16_xz(:,:,itMax) = cylToCart(Ipropagate16(:,:,itMax),r,x);
% [~, Ipropagate64(:,:,itMax)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex64(:,:,itMax));
% Ipropagate64_xz(:,:,itMax) = cylToCart(Ipropagate64(:,:,itMax),r,x);
% [~, Ipropagate128(:,:,itMax)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex128(:,:,itMax));
% Ipropagate128_xz(:,:,itMax) = cylToCart(Ipropagate128(:,:,itMax),r,x);
% 
% 
% PF_probeComparison(Iprobe_xz, Ipropagate16_xz(:,:,itMax), x, z, t(itMax), izMax);
% figure;
% hold on;
% plot(x*1e6, Iprobe_xz(:,1), "r", LineWidth=2, LineStyle= ":");
% plot(x*1e6, Ipropagate16_xz(:,102,itMax), "m", LineWidth=2);
% plot(x*1e6, Ipropagate_xz(:,izMax,itMax), "y", LineWidth=2, LineStyle= '-.');
% plot(x*1e6, Ipropagate64_xz(:,56,itMax), "c", LineWidth=2, LineStyle= '--');
% plot(x*1e6, Ipropagate128_xz(:,39,itMax), "b", LineWidth=2,  LineStyle= '-');
% hold off;
% axis tight; grid on;
% xlabel('x [\mum]');ylabel('Intensity [W/m^2]');
% legend('0[nJ]','16[nJ]','35[nJ]','64[nJ]','128[nJ]');
% title('Maximal Probe Intensity with Different Pump Energies');

% Change probe's w0:


% Change probe's pulse duration:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

