%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN SIMULATION FILE
% -------------------------------------------------------------------------
% Defines all the necessary structures, conducts the experiment flow, plots
% graphs and calculates necessary tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paraneters and Vector Definitions
% =========================================================================

sp = systemParameters();

% Temporal vector intialization:
tf = 3000;                              % 3[ns]
t = 0:10:tf;
t = t .* 1e-12;                         % Time vector, 0-3[ns], [s]
Nt = numel(t);

% Spatial vectors intialization:
Nr = 501;                               % Number of elements
Nz = 501;

z = linspace(0, sp.Lz2, Nz);            % 25 micron

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
pump_width = 3e-11;     % Pulse of 30[ps]
pump_E = 35e-9;         % Beam total energy [J]
pump_w0 = 2.7e-6;                       
z0 = 5e-6;                 % Beam waist location

% Pump laser:
pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, pump_w0, z0);  % 775[nm]
Ipump = pump.intensityProfileBLDumped(z);
Ipump_xz = cylToCart(Ipump,r,x);

% Plotting pump intensity:
PF_plot_xz(Ipump_xz*1e-4, z, x, "Pump Intensity [W/cm^2]");
%PF_plot_xz(Ipump_xz/max(Ipump_xz(:)), z, x, "Pump Normalized");
%PF_x(Ipump_xz(:,1)*1e-4,x,"Pump Intensity at Surface (z=0)","Intensity [W/cm^2]",'Intensity [W/cm^2]')

iz = 1;     % z=0
%PF_beamIntensityPlot_xy(Ipump(:,iz), r, x, z(iz));
%PF_beamIntensityPlot_xy(Ipump(:,1)/max(Ipump(:,1)), r, x, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FCC Creation & Diffusion
% =========================================================================

% Shining pump beam on the sample, causing e-h generation and diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)

% Finding the time when we get maximal concentration:
[~, ~, itMax] = findMax(pDiff);   % itMax = 12 = 110[ps]


pxz = cylToCart(pDiff(:,:,itMax),r,x) * 1e-6;
%PF_x(pxz,x,"FCC Maximal Concentration in [1/cm^3] at t=110[ps]")
%PF_plot_xz(pxz,z,x,"FCC Maximal Concentration in [1/cm^3] at t=110[ps]");


%tIdx = [1, 19, 20, 50, 100, 150, 270, Nt];
%PF_colormapAnimation_xz(pDiff*1e-6, r, z ,x, 1:Nt, "FCC Concentration vs. Time", "FCC concentration [1/cm^3]", false, 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refractive Index Calculation
% =========================================================================

% Calculating the complex refractive index n(r,z,t) (changes due to FCC generation):
n_complex = complexRefractiveIndex(pDiff, pump.lambda);

% Finding the z index where we get maximal absorption:
[~, izMax, ~] = findMax(imag(n_complex(:,:,itMax)));

%PF_complexRefractiveIndex(n_complex(:,:,itMax), r, z, x, pump.lambda, izMax, t(itMax));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probe Beam
% =========================================================================

% Probe laser:
probe_wd = 5e-11;   % Pulse of 50[ps]
probe_E = pump.pulse_energy/100;        % [J]
probe_w0 = 2.7e-6;                       

probe = laser(sp.wl2, probe_wd, probe_E, "Gauss", r, phi, z, probe_w0, 0);   % 775[nm]
Iprobe = probe.intensityProfileBLDumped(z);
Iprobe_xz = cylToCart(Iprobe,r,x);

%PF_plot_xz(Iprobe_xz, z, x, "Probe Intensity [W/m^2]");
%PF_plot_xz(Iprobe_xz/max(Iprobe_xz(:)), z, x, "Probe Normalized");
PF_x(Iprobe_xz(:,1),x,"Undisturbed Probe Intensity at Surface (z=0)");

%absorptionDepth(Iprobe(1,:),z,probe.lambda);

% Cutting the Probe:
%aperture = rectpuls(r,3.506e-6);
%Icut = aperture' .* Iprobe(:,1);
%Icut_xz = cylToCart(Icut,r,x);
%PF_x(Icut_xz,x,"Probe Intensity at Surface After Cutting");

%PF_FWHM_z(FWHM(Iprobe,r),z,"Undisturbed Probe FWHM");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FCC Slicing
% =========================================================================

% Creating FCC slices (FCC only at specific z, ni elsewhere):
%slice_width = 0.5;  % In microns. 25 = 1 slice

%[deltaFWHM_all, minFWHM_all, izSlices] = sliceContribution(pDiff(:,:,itMax), pump, probe, r, z, slice_width);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probe Propogation (BPM)
% =========================================================================

% Ipropogate = zeros(Nr,Nz,Nt);       % Probe intensity after propogation
% Ipropogate_xz = zeros(Nx,Nz,Nt);
% 
% for i = 1:Nt
%     [~, Ipropogate(:,:,i)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,i));
%     %Ipropogate_xz(:,:,i) = cylToCart(Ipropogate(:,:,i),r,x);
% end

% No time coordinate:
[~, I] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,itMax));
I_xz = cylToCart(I,r,x);

% Finding the z index where we get maximal intensity or focusing:
% [~, izMax,~] = findMax(Ipropogate_xz(:,:,itMax));
[~, izMax,~] = findMax(I_xz);
prpFWHM = FWHM(I(:,:),r);
[~,i] = min(prpFWHM);

%PF_plot_xz(Ipropogate_xz(:,:,itMax),z,x)
%PF_x(Ipropogate_xz(:,izMax,itMax),x,"Probe Intensity at Maximum Focusing and 110[ps] delay");
%PF_plot_xz(I_xz,z,x,"Probe Intensity 110[ps] delay",izMax);
% PF_x_alongz(Iprobe_xz,x,z,"Undisturbed Probe Intensity 110[ps] delay");
PF_x(I_xz(:,i),x,"Probe Intensity at Maximum Focusing and 110[ps] delay");
PF_x_alongz(I_xz,x,z,"Probe Intensity 110[ps] delay");
PF_x(I_xz(:,end),x,"Probe Intensity at sample end and 110[ps] delay");

%PF_FWHM_z(FWHM(Iprobe,r),z,"Undisturbed Probe FWHM");
% PF_FWHM_z(prpFWHM,z,"Propagated Probe FWHM");

%nxz = cylToCart(n_complex(:,:,itMax),r,x);
%PF_x(imag(nxz(:,:)),x,"k 110[ps] delay");
%PF_x(real(nxz(:,:)),x,"n 110[ps] delay");

%PF_colormapAnimation_xz(Ipropogate*1e-4,r,z,x,1:Nt,"Probe Propogation", "Intensity [W/cm^2]", false, 50)

%PF_x(Ipropagate_xz(201,:,itMax),z,"Probe Intensity along x=0 at 110[ps] delay");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons and Other Tests
% =========================================================================

% Checking the best z0 and w0:
z0_vec = (6:1:10)*1e-6;
w0_vec = (1.17:0.005:1.2)*1e-6;
%CheckZ0(pump,z,r,x,t,probe);

% Checking the best pump energy:
%[opt_energy,opt_fwhm] = optimalEnergy(pump,z,r,x,t,probe);

% FWHM calculation:
%fwhm = zeros(Nz,Nt);

%for i = 1:Nt
%    fwhm(:,i) = FWHM(Ipropogate(:,:,i),r);
    %fwhm(:,i) = FWHM(Ipropagate_xz(:,:,i),x);
%end

%PF_FWHM_z(fwhm(:,itMax),z,sprintf('Propogated Probe FWHM at %d[ps]', t(itMax)*1e12));

% FWHM vs. t:
%tmin = min(fwhm(:));
% figure;
% plot(t,tmin(:))
% grid on; axis tight;
% xlabel("t [ps]"); ylabel('Minimal FWHM [\mum]');

% Comparing the probe with and without pump:I = {Iprobe_xz, Ipropogate_xz(:,:,itMax)};      % List of intensities Nx x Nz
names = {'0[nJ]','35[nJ]'};
izList = {i,i};
PF_probeComparison({Iprobe_xz,I_xz}, names, x, z, t(itMax), izList);

% Checking the maximal intensity vs. delay time:
%PF_maxIntensity(Ipropagate_xz, t);   


%PF_x(Ipropagate_xz(201,:,itMax),z,"Probe Intensity along x=0 at itmax")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pump Energy Related Effects
% =========================================================================

% Different pump energies:
% energies = [0, 10, 20, 30, 35, 64] * 1e-9;   % [J]
% names = {'0[nJ]','10[nJ]','20[nJ]','30[nJ]','35[nJ]','64[nJ]'};
% 
% Ilist = cell(1, numel(energies));
% izList = cell(1, numel(energies));
% 
% for k = 1:numel(energies)
% 
%     if energies(k) == 0
%         % No pump:
%         probeK = laser(probe.lambda, probe.pulse_width, probe.pulse_energy, "Gauss", r, phi, z, probe.w0, probe.z0);
% 
%         [~, I_rz] = propagationBPM_rz(probeK.profile(:,1), r, z, probeK, ones(Nr,Nz)*(sp.n+1i*(sp.alpha*probe.lambda/(4*pi))));
% 
%     else
%         pumpK = laser(sp.wl2, pump_width, energies(k), "Donut", r, phi, z, w0, z0);
% 
%         probeK = laser(probe.lambda, probe.pulse_width, energies(k)/100, "Gauss", r, phi, z, probe.w0, probe.z0);
% 
%         n_complexK = complexRefractiveIndex(FCCDiffusion(pumpK, t, r, z), pumpK.lambda);
% 
%         [~, I_rz] = propagationBPM_rz(probeK.profile(:,1), r, z, probeK, n_complexK(:,:,itMax));
%     end
% 
%     I_xz = cylToCart(I_rz, r, x);
%     PF_x_alongz(I_xz,x,z,names(k));
% 
%     %[~, izK] = findMax(I_xz);
%     %PF_x(I_xz(:,izK),x,names(k)+" max intensity");
% 
%     fwhm_z = FWHM(I_rz,r);
%     PF_FWHM_z(fwhm_z,z,names(k));
% 
%     [~,izK] = min(fwhm_z);
%     PF_x(I_xz(:,izK),x,names(k)+" min fwhm");
%     PF_x(I_xz(:,end),x,names(k)+" end");
% 
%     Ilist{k} = I_xz;
%     izList{k} = izK;
% end
% 
% PF_probeComparison(Ilist, names, x, z, t(itMax), izList);

% Change probe's w0:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

