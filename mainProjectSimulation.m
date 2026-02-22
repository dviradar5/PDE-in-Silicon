%% DOCUMENT
%% CLEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Simulation file
% -------------------------------------------------------------------------
% defines all the necessary structures, conducts the experiment flow, plots
% graphs and calculates necessary 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paraneters and vector definitions:
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

r = linspace(0,sqrt(2)*sp.Lx,Nr);       % 10x10 micron sample so r covers it
phi = atan(1);                          % y = x, radial symmetry

x = linspace(-max(r), max(r), 2*Nr-1);
Nx = numel(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ENERGIES, INTENSITY OoM, beam diameter
% Pump parameters:
pmp_wd = 3e-11;         % Pulse of 30[ps]
pmp_D = sp.D2;
w0 = 5e-6;
z0 = 0;                 % Beam waist location
pump_E = 35e-9;         % Beam total energy [J]       how does it get into account????????

% Pump laser:
pump = laser(sp.wl2, pmp_wd, pump_E, pmp_D, "Donut", r, phi, z, w0, z0);  % 775[nm]
Ipump = pump.intensityProfileBLDumped(z);
Ipump_xz = cylToCart(Ipump,r,x);

% Plotting pump intensity:
%PF_beamIntensityPlot_xz(Ipump_xz, z, x, "Pump Intensity [W/m^2]");
%PF_beamIntensityPlot_xz(Ipump_xz/max(Ipump_xz(:)), z, x, "Pump Normalized");

iz = 1;     % z=0
%PF_beamIntensityPlot_xy(Ipump(:,iz), r, x, z(iz));
%PF_beamIntensityPlot_xy(Ipump(:,1)/max(Ipump(:,1)), r, x, 1);

% NOTE: pump doesn't need to be t-dependant. Time in this simulation can be
% regarded as the time since when the pump hit the sample.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shining pump beam on the sample, causing e-h generation and diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)

tIdx = [1, 19, 20, 50, 100, 150, 270, Nt];
%PF_colormapAnimation_xz(pDiff*1e-6, r, z ,x, 1:Nt, "FCC Concentration vs. Time", "FCC concentration [1/cm^3]", false, 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% absorption OoM
% Calculating the complex refractive index n(r,z,t) (changes due to FCC generation):
n_complex = complexRefractiveIndex(pDiff, pump.lambda);

it = 12;
PF_complexRefractiveIndex(n_complex(:,:,it), r, z, x, pump.lambda, iz, t(it));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ENERGY, GAUSSIAN PULSE EVELOPE
% Probe laser:
probe_wd = 5e-11;   % pulse of 50[ps]  how does it affect calculation?????
probe_D = sp.D1;
probe_E = pump.pulse_energy/100;         % [J] 
probe_w0 = w0;                          % how does it affect calculation?????

probe = laser(sp.wl2, probe_wd, probe_E, probe_D, "Gauss", r, phi, z, probe_w0, 0);   % 775[nm]
Iprobe = probe.intensityProfileBLDumped(z);
Iprobe_xz = cylToCart(Iprobe,r,x);

PF_beamIntensityPlot_xz(Iprobe_xz, z, x, "Probe Intensity [W/m^2]");
%PF_beamIntensityPlot_xz(Iprobe_xz/max(Iprobe_xz(:)), z, x, "Probe Normalized");

%PF_FWHM_z(FWHM(Iprobe,r),z,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probe propogation (BPM):
Ipropagate = zeros(Nr,Nz,Nt);           % Probe intensity after propogation
Ipropagate_xz = zeros(Nx,Nz,Nt);

for i = 1:Nt
    [~, Ipropagate(:,:,i)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,i));
    Ipropagate_xz(:,:,i) = cylToCart(Ipropagate(:,:,i),r,x);
end

%PF_beamIntensityPlot_xz(Ipropagate_xz(:,:,11), z, x, "Probe Intensity with 0 delay");
%PF_beamIntensityPlot_xz(Ipropagate_xz/max(Ipropagate_xz(:)), z, x, "Probe Normalized");
%PF_intensity_x(Ipropagate_xz(:,120,11),x);
PF_colormapAnimation_xz(Ipropagate,r,z,x,1:Nt,"Probe Propogation", "Intensity [W/m^2]", false, 50)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparisons and other calculations:

% Comparing the probe with and without pump:
it = 11;
%PF_probeComparison(Iprobe_xz, Ipropagate_xz(:,:,it), x, z, t(it));

% Checking the maximal intensity vs. delay time:
%PF_maxIntensity(Ipropagate_xz, t);   

% FWHM calculation:
fwhm = zeros(Nz,Nt);

for i = 1:Nt
    fwhm(:,i) = FWHM(Ipropagate(:,:,i),r);
end

it = 11;
PF_FWHM_z(fwhm(:,it),z,t(it));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pump energy related effects:

