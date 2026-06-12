%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation for Asaf's Paper
% -------------------------------------------------------------------------
% In this file we run the simulation with the parameters from Asaf's paper,
% with the goal of seeing good correlation between the results
% 
% We assume the Beam Expender was deployed and use the parameter values
% that it created, this way we don't need to implement it into the
% simulation, only check whether the BE-shrunk pump values corresponds to
% the probe values
%
% Maybe its wrong to think that way, since the beam expender also made it
% denser and more parallel, so the effects along the sample are also
% different
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paraneters and Vector Definitions
% =========================================================================

sp = systemParameters();

% Temporal vector intialization:
tf = 1000;                              % 1[ns] in [ps]
t = 0:10:tf;
t = t .* 1e-12;                         % Time vector, 0-1[ns], [s]
Nt = numel(t);

% Spatial vectors intialization:
Nr = 501;                               % Number of elements
Nz = 501;

z = linspace(0, sp.Lz1, Nz);            % 20 micron

r = linspace(0,sp.Lx,Nr);               % 10x10 micron sample
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
pump_width = 30e-12;                    % Pulse of 30[ps]
pump_E = 20e-9;                         % Beam total energy of 20[nJ], [J]
pump_w0 = 1.185e-6;%1.185e-6; %3.345e-6
z0 = 9e-6;%20e-6                        % Beam waist location
l = 1;                                  % LG polynomial index

% Pump laser:
pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, pump_w0, z0,l);  % 775[nm]

Ipump = pump.intensityProfileBLDumped(z);
Ipump_xz = cylToCart(Ipump,r,x);

% Plots:
[~,f1] = PF_x(Ipump_xz(:,end)*1e-4,x,"Pump Intensity at Sample End",'Intensity [W/cm^2]');
%[~,f1] = PF_x(Ipump_xz(:,end)/max(Ipump_xz(:,end)),x,'','Intensity');
%PF_FWHM_z(FWHM(Ipump,r),z)                                         %FIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probe Beam
% =========================================================================

% Probe parameters:
probe_width = 30e-12;                   % Pulse of 30[ps]
probe_E = pump.pulse_energy/100;        % [J]
probe_w0 = 2.28e-6;                       

% Probe laser:
probe = laser(sp.wl2, probe_width, probe_E, "Gauss", r, phi, z, probe_w0, 0);   % 775[nm]

Iprobe = probe.intensityProfileBLDumped(z);
Iprobe_xz = cylToCart(Iprobe,r,x);

% Plots:
[~,f2] = PF_x(Iprobe_xz(:,end)*1e-4,x,"Undisturbed Probe Intensity at the Sample End",'Intensity [W/cm^2]');
%[~,f2] = PF_x(Iprobe_xz(:,end)/max(Iprobe_xz(:,end)),x,'','Intensity');
% PF_FWHM_z(FWHM(Iprobe,r),z,"Undisturbed Probe FWHM");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters Optimization and Other Tests
% =========================================================================

% Checking for the optimal waist, waist location and energy:
% z0_vec = (0:1:20)*1e-6;%(-10:5:30)*1e-6;
% w0_vec = (1.1:0.1:1.4)*1e-6;

% CheckZ0(pump,z,r,x,t,probe,z0_vec,l);
% CheckW0andZ0(pump,z,r,x,t,probe,z0_vec,w0_vec,0.91e-6,0.96e-6);

% [opt_energy,opt_fwhm] = optimalEnergy(pump,z,r,x,t,probe);

% Checking different pump vortex orders:
l_vec = 1:1:8;
%CheckLGOrder(pump,z,r,x,t,probe,l_vec);

% Checking we're under the destructive intensity:
damageThreshold(Ipump, pump_width, pump.lambda, 'Pump:');
damageThreshold(Iprobe, probe_width, probe.lambda, 'Undisturbed Probe:');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Run
% =========================================================================

% FCC diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)
[~, ~, itMax] = findMax(pDiff);

% Calculating the new complex refractive index:
n_complex = complexRefractiveIndex(pDiff, pump.lambda);
[~, in] = findMax(n_complex(:,:,itMax));% Maximum intensity

% PML mask parameters:
a = 4;
b = 8;

%[pDiff,n_complex,Ipropagate,Ipropagate_xz] = runSimulation(x,z,r,t,pump,probe,a,b);

% Simulation for maximal delay (~90-110[ps]):
[~, Ipropagate] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,itMax),a,b);
Ipropagate_xz = cylToCart(Ipropagate,r,x);

% Checking we're under the destructive intensity:
damageThreshold(Ipropagate, probe_width, probe.lambda, 'Propagated Probe:');

[~,f3] = PF_x(Ipropagate_xz(:,Nz)*1e-4,x,"Propagated Probe Intensity at the Sample End",'Intensity [W/cm^2]');
%[~,f3] = PF_x(Ipropagate_xz(:,end)/max(Ipropagate_xz(:,end)),x,'','Intensity');
%PF_x_alongz(Ipropagate_xz*1e-4, x, z,"Propagated Probe Intensity",'Intensity [W/cm^2]');

% Combining the figures:
% figure('Color','w');
% ax1 = subplot(1,3,1);
% copyobj(allchild(findobj(f2,'Type','axes')), ax1);
% ylabel('Normalized Intensity', 'FontSize',12); xlabel('x [\mum]');
% 
% ax2 = subplot(1,3,2);
% copyobj(allchild(findobj(f1,'Type','axes')), ax2);
% xlabel('x [\mum]');
% 
% ax3 = subplot(1,3,3);
% copyobj(allchild(findobj(f3,'Type','axes')), ax3);
% ylim([0,1]);
% xlabel('x [\mum]');

[~, i] = findMax(Ipropagate_xz);        % Maximum intensity

prpFWHM = FWHM(Ipropagate,r);
%[~,i] = min(prpFWHM);                  % Minimal FWHM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other Plots
% =========================================================================

% Pump:
PF_x_alongz(Ipump_xz*1e-4, x, z,"Pump Intensity",'Intensity [W/cm^2]');
PF_x(Ipump_xz(:,1)*1e-4,x,"Pump Intensity at the Sample Surface (z=0)", "Intensity [W/cm^2]");
PF_beamIntensityPlot_xy(Ipump(:,1)*1e-4, r, x);

% Original probe:
PF_x_alongz(Iprobe_xz*1e-4, x, z,"Undisturbed Probe Intensity",'Intensity [W/cm^2]');
PF_x(Iprobe_xz(:,1)*1e-4,x,"Undisturbed Probe Intensity at the Sample Surface (z=0)",'Intensity [W/cm^2]');
PF_x(Iprobe_xz(:,i)*1e-4,x,"Undisturbed Probe Intensity at Maximum Intensity",'Intensity [W/cm^2]',[]);
PF_beamIntensityPlot_xy(Iprobe(:,1)*1e-4, r, x,"Undisturbed Probe Intensity at the Sample Surface (z=0) [W/cm^2]");

% Complex refractive index:
PF_complexRefractiveIndex(n_complex(:,:,itMax), r, z, x, pump.lambda, 1, t(itMax));

% Propagated probe:
PF_x_alongz(Ipropagate_xz*1e-4, x, z,"Propagated Probe Intensity",'Intensity [W/cm^2]');
PF_x(Ipropagate_xz(:,i)*1e-4,x,"Propagated Probe Intensity at Maximum Intensity",'Intensity [W/cm^2]');
PF_x(Ipropagate_xz(:,Nz)*1e-4,x,"Propagated Probe Intensity at the Sample End",'Intensity [W/cm^2]');
PF_beamIntensityPlot_xy(Ipropagate(:,i)*1e-4, r, x, 'Propagated Probe Intensity at Maximum Intensity [W/cm^2]');

% FWHM:
PF_FWHM_z(prpFWHM,z,"Propagated Probe FWHM");
