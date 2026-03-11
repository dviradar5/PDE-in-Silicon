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
[~, ~, itMax] = findMax(pDiff);   % itMax = 12 = 110[ps]

%pxz = cylToCart(pDiff(:,:,itMax),r,x) * 1e-6;
%PF_plot_xz(pxz,z,x,"FCC Maximal Concentration in [1/cm^3] at t=110[ps]");


%tIdx = [1, 19, 20, 50, 100, 150, 270, Nt];
%PF_colormapAnimation_xz(pDiff*1e-6, r, z ,x, 1:Nt, "FCC Concentration vs. Time", "FCC concentration [1/cm^3]", false, 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refractive Index Calculation
% -------------------------------------------------------------------------

% Calculating the complex refractive index n(r,z,t) (changes due to FCC generation):
n_complex = complexRefractiveIndex(pDiff, pump.lambda);

% Finding the z index where we get maximal absorption:
[~, izMax, ~] = findMax(imag(n_complex(:,:,itMax)));

%PF_complexRefractiveIndex(n_complex(:,:,itMax), r, z, x, pump.lambda, izMax, t(itMax));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probe Propogation (BPM)
% -------------------------------------------------------------------------

Ipropagate = zeros(Nr,Nz,Nt);       % Probe intensity after propogation
Ipropagate_xz = zeros(Nx,Nz,Nt);

for i = 1:Nt
    [~, Ipropagate(:,:,i)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,i));
    Ipropagate_xz(:,:,i) = cylToCart(Ipropagate(:,:,i),r,x);
end

% Finding the z index where we get maximal intensity:
[~, izMax,~] = findMax(Ipropagate_xz(:,:,itMax));

%PF_intensity_x(Ipropagate_xz(:,izMax,itMax),x,"Probe Intensity at Maximum Focusing and 110[ps] delay");
%PF_colormapAnimation_xz(Ipropagate,r,z,x,1:Nt,"Probe Propogation", "Intensity [W/m^2]", false, 50)

%PF_intensity_x(Ipropagate_xz(201,:,itMax),z,"Probe Intensity along x=0 at 110[ps] delay");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons and Other Calculations
% -------------------------------------------------------------------------

% Comparing the probe with and without pump:
I = {Iprobe_xz, Ipropagate_xz(:,:,itMax)};      % List of intensities NxxNz
names = {'0[nJ]','35[nJ]'};
izList = {izMax,izMax};
%PF_probeComparison(I, names, x, z, t(itMax), izList);

% Checking the maximal intensity vs. delay time:
%PF_maxIntensity(Ipropagate_xz, t);   

% FWHM calculation:
fwhm = zeros(Nz,Nt);

for i = 1:Nt
    fwhm(:,i) = FWHM(Ipropagate(:,:,i),r);
    %fwhm(:,i) = FWHM(Ipropagate_xz(:,:,i),x);
end

PF_FWHM_z(fwhm(:,itMax),z,t(itMax));

%PF_intensity_x(Ipropagate_xz(201,:,itMax),z,"Probe Intensity along x=0 at itmax")

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
% I16 = zeros(Nr,Nz,Nt); I64 = zeros(Nr,Nz,Nt); I128 = zeros(Nr,Nz,Nt);
% I16_xz = zeros(Nx,Nz,Nt); I64_xz = zeros(Nx,Nz,Nt); I128_xz = zeros(Nx,Nz,Nt);
% 
% [~, I16(:,:,itMax)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex16(:,:,itMax));
% I16_xz(:,:,itMax) = cylToCart(I16(:,:,itMax),r,x);
% [~, I64(:,:,itMax)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex64(:,:,itMax));
% I64_xz(:,:,itMax) = cylToCart(I64(:,:,itMax),r,x);
% [~, I128(:,:,itMax)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex128(:,:,itMax));
% I128_xz(:,:,itMax) = cylToCart(I128(:,:,itMax),r,x);
% 
% [~, iz16,~] = findMax(I16_xz(:,:,itMax));
% [~, iz64,~] = findMax(I64_xz(:,:,itMax));
% [~, iz128,~] = findMax(I128_xz(:,:,itMax));
% 
% I = {Iprobe_xz, I16_xz(:,:,itMax), Ipropagate_xz(:,:,itMax), I64_xz(:,:,itMax), I128_xz(:,:,itMax)};
% names = {'0[nJ]','16[nJ]','35[nJ]','64[nJ]','128[nJ]'};
% izList = {1,iz16,izFocus,iz64,iz128};
% PF_probeComparison(I, names, x, z, t(itMax), izList);
% 
% PF_intensity_x(Iprobe_xz(201,:),z,"Probe Intensity along x=0  0")
% PF_intensity_x(I16_xz(201,:,itMax),z,"Probe Intensity along x=0 at 110[ps] delay 16")
% PF_intensity_x(I64_xz(201,:,itMax),z,"Probe Intensity along x=0 at 110[ps] delay 64")
% PF_intensity_x(I128_xz(201,:,itMax),z,"Probe Intensity along x=0 at 110[ps] delay 128")

% Change probe's w0:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

