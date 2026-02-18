%% DOCUMENT
clc; clear;
%close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paraneters and vector definitions:
sp = systemParameters();

% Temporal vector intialization:
tf = 30000;                 % 30[ns]
t = 0:10:tf;
t = t .* 1e-12;             % Time vector, 0-30[ns], [s]
Nt = numel(t);

% Spacial vector intialization:
Nr = 201;                               % Number of elements
Nz = 201;

z = linspace(0, sp.Lz, Nz);

r = linspace(0,sqrt(2)*sp.Lx,Nr);       % 10x10 micron sample so r covers it
phi = atan(1);                          % y = x, radial symmetry

x = linspace(-max(r), max(r), 2*Nr-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ENERGIES, INTENSITY OoM
% Laser parameters:
pmp_wd = 3e-11;         % Pulse of 30[ps]
pmp_D = sp.D2;
w0 = 5e-6;
z0 = 0;                 % Beam waist location (-7 in old simulation)
pump_E = 35e-9;         % [J] FIXXXXXXXXXXXXX

% Pump laser:
pump = laser(sp.wl2, pmp_wd, pump_E, pmp_D, "Donut", r, phi, z, w0, z0);  % 775[nm]

Ipump = pump.intensityProfileBLDumped(z);

PF_beamIntensityPlotXZ(Ipump, r, z, x)
%PF_beamPlotXY(r, x, Ipump(:,1), 1, 10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shining pump beam on the sample, causing e-h generation and diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution
tIdx = [1, 20, 200, 1000, Nt];
%PF_colormapAnimation_xz(pDiff*1e-6, r, z ,x, tIdx, "FCC Concentration vs. Time", "FCC concentration [1/cm^3]");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAYBE MOVE TO FUNCTION OR MERGE WITH EPSILON
% Calculating the complex refractive index n(r,z,t) (changes due to FCC generation):
n_complex = zeros(Nr, Nz, Nt);      % complex refractive index over time

for it = 1:Nt
    N = pDiff(:,:,it);              % Ne = Nh since 1 photon = e+h

    % Bennett–Soref expects Ne, Nh at pump wavelength:
    [dn, dalpha] = BennettSoref(N, N, pump.lambda);

    n_complex(:,:,it) = (sp.n + dn) + 1i*(sp.alpha + dalpha)*pump.lambda/(4*pi);
    
    % Calculating relative permitivity for maxwell's equations later:
    %epsilon_r = n_complex^2;
end

PF_complexRefractiveIndex(n_complex(:,:,1), r, z, x, pump.lambda, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ENERGY, GAUSSIAN PULSE EVELOPE
% Probe laser:
prb_wd = 5e-11;   % pulse of 50[ps]
prb_D = sp.D1;
probe_E = pump.pulse_energy/100;         % [J]

probe = laser(sp.wl2, prb_wd, probe_E, prb_D, "Gauss", r, phi, z, 3*w0, 0);   % 775[nm]

Iprobe = probe.intensityProfileBLDumped(z);

PF_beamIntensityPlotXZ(Iprobe,r,z,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probe propogation (BPM):
Iproagate = zeros(Nr,Nz,numel(tIdx));           % Probe intensity after propogation

for i = 1:numel(tIdx)
    [~, Iproagate(:,:,i)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,i));

    %PF_beamIntensityPlotXZ(I_rz, r, z, x);
end

%PF_colormapAnimation_xz(Iproagate,r,z,x,tIdx,"Probe Propogation", "Intensity [W/m^2]")

% Comparing the probe with and without pump:
PF_probeComparison();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


