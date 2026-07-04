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
tf = 150;                               % 150 [ps]
t = 0:5:tf;
t = t .* 1e-12;                         % Time vector, 0-300 [ps], [s]
Nt = numel(t);

% Spatial vectors intialization:
Nr = 1001;                              % Number of elements
Nz = 1001;

z = linspace(0, sp.Lz2, Nz);            % 25 micron

r = linspace(0,4*sp.Lx,Nr);             % 20x20 micron sample
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
pump_width = 30e-12;                    % Pulse of 30[ps] 50
pump_E = 35e-9;                         % Beam total energy [J]
pump_w0 = 4.24e-6;                      % 1.185 3.385
z0 = 25e-6;                             % Beam waist location
l = 1;                                  % LG polynomial index

% Probe parameters:
% probe_width = 30e-12;                   % Pulse of 50[ps]
% probe_w0 = 5e-6;                    % 1.875e-6          try really wide one           
% n = (sp.n + 1i*sp.alpha*sp.wl2/(4*pi))*ones(Nr,Nz);                     % Includes BL 
% 
% probe = laser(sp.wl2, probe_width, pump_E/100, "Gauss", r, phi, z, probe_w0, 0, sp.n);   % 775[nm]
% 
% Different energies:
% L = [0,50e-9,75e-9,100e-9];
% labels = ["No Pump","50 [nJ]","75 [nJ]","100 [nJ]"];
% L = [0,1.185e-6,3.385e-6,5e-6,7.5e-6];%,1.185e-6
% labels = ["No Pump","1.185 [\mum]","3.385 [\mum]","5 [\mum]","7.5 [\mum]"];%,"1.185 [\mum]"
% L = [1,2];
% labels = ['1','2'];
% 
% Dfigs = gobjects(1,numel(L));
% Ifigs = gobjects(1,numel(L));
% 
% for i = 1:numel(L)
%     %pump = laser(sp.wl2, pump_width, L(i), "Donut", r, phi, z, pump_w0, z0, sp.n, l);
%     %pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, L(i), z0, sp.n, l);
%     pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, pump_w0, z0, sp.n, L(i));
% 
%     %Ipump = pump.intensityProfileBLDumped(z);
%     %Ipump_xz = cylToCart(Ipump,r,x);
%     %PF_plot_xz(Ipump_xz/max(Ipump_xz(:)), z, x, "Pump Normalized");
% 
%     pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)
%     [~, ~, itMax] = findMax(pDiff); 
% 
%     n_complex = complexRefractiveIndex(pDiff, pump.lambda);
%     %PF_complexRefractiveIndex(n_complex(:,:,itMax), r, z, x, pump.lambda, 1, t(itMax));
% 
%     if L(i) == 0
%         %probe = laser(sp.wl2, probe_width, 1e-9, "Gauss", r, phi, z, probe_w0, 0, sp.n);
%         [~, I] = propagationBPM_rz(probe.profile(:,1),r,z,probe,n,0,8);    % [W/m^2]
%     else
%         %probe = laser(sp.wl2, probe_width, L(i)/100, "Gauss", r, phi, z, probe_w0, 0, sp.n);
%         [~, I] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,itMax),4,20);
%     end
% 
%     I_xz = cylToCart(I,r,x);
% 
%     %interactiveIxMovie(I_xz*1e-4, x, z, "Probe I(x,z)");
%     %PF_plot_xz(I_xz*1e-4, z, x, "Probe Intensity");
% 
%     Imax = max(I_xz);% .* exp(sp.alpha * z);
%     Ifigs(i) = figure('Color','w');
%     plot(z*1e6,Imax/Imax(1),'LineWidth', 3);
%     xlabel("z [\mum]"); ylabel("Intensity [au]");
%     title(labels(i));
% 
%     %PF_x(I_xz(:,end)*1e-4,x,"Undisturbed Probe Intensity at the Sample End","Intensity [W/cm^2]");
%     PF_x_alongz(I_xz*1e-4,x,z,"Propagated Probe Intensity");
% 
%     [Dfigs(i),~] = computeDiameter(I, r, z, "Gauss","Propogated Probe Diameter"); 
% end
% 
% PF_combinePlots(L, labels, Dfigs,"Width [\mum]" , '');
% PF_combinePlots(L, labels, Ifigs,"Intensity [au]" , '');

% Pump laser:
pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, pump_w0, z0, sp.n, l);  % 775[nm]

Ipump = pump.intensityProfileBLDumped(z);
Ipump_xz = cylToCart(Ipump,r,x);

% Plotting pump intensity:
%PF_plot_xz(Ipump_xz/max(Ipump_xz(:)), z, x, "Pump Normalized");
%PF_x(Ipump_xz(:,end)*1e-4,x,"Pump Intensity at the Sample End","Intensity [W/cm^2]");
%PF_x(Ipump_xz(:,1)*1e-4,x,"Pump Intensity at the Sample Front","Intensity [W/cm^2]");

%computeDiameter(Ipump, r, z, "Donut","Pump Diameters");
damageThreshold(Ipump, pump_width, pump.lambda, "Pump:");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FCC Creation & Diffusion
% =========================================================================

% Shining pump beam on the sample, causing e-h generation and diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)

% Finding the time when we get maximal concentration:
[~, ~, itMax] = findMax(pDiff);        % itMax = 12 = 110[ps]
%pxz = cylToCart(pDiff(:,:,itMax),r,x) * 1e-6;

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
probe_width = 30e-12;                   % Pulse of 30[ps]
probe_E = pump.pulse_energy/100;        % [J]
probe_w0 = 2.267e-6;                    % 2.378 [um]
z0 = 0;                                 % sp.Lz2

probe = laser(sp.wl2, probe_width, probe_E, "Gauss", r, phi, z, probe_w0, z0, sp.n);   % 775[nm]

% Undisturbed probe propagating through the sample:
n = (sp.n + 1i*sp.alpha*sp.wl2/(4*pi))*ones(Nr,Nz);                     % Includes BL 
[~, Iprobe] = propagationBPM_rz(probe.profile(:,1),r,z,probe,n,0,8);    % [W/m^2]
%Iprobe = probe.intensityProfileBLDumped(z);
Iprobe_xz = cylToCart(Iprobe,r,x);

%PF_plot_xz(Iprobe_xz/max(Iprobe_xz(:)), z, x, "Probe Normalized");
%PF_x(Iprobe_xz(:,1)*1e-4,x,"Undisturbed Probe Intensity at the Sample Front","Intensity [W/cm^2]");
PF_x(Iprobe_xz(:,end)*1e-4,x,"Undisturbed Probe Intensity at the Sample End","Intensity [W/cm^2]");
PF_x_alongz(Iprobe_xz*1e-4,x,z,"Propagated Probe Intensity");

%absorptionDepth(Iprobe(1,:),z,probe.lambda);
damageThreshold(Iprobe, probe_width, probe.lambda, "Undisturbed Probe:");

% Probe Diameter:
%computeDiameterTheoretically(z, "Gauss",0, sp.wl2, probe_w0, 0, sp.n);
%computeDiameter(Iprobe, r, z, "Gauss");
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

% PML mask parameters:
a = 4;
b = 8;

% Ipropogate = zeros(Nr,Nz,Nt);          % Probe intensity after propogation
% Ipropogate_xz = zeros(Nx,Nz,Nt);
% 
% for i = 1:Nt
%     [~, Ipropogate(:,:,i)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,i),a,b);
%     %Ipropogate_xz(:,:,i) = cylToCart(Ipropogate(:,:,i),r,x);
% end

% No time coordinate:
[~, I] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,itMax),a,b);
I_xz = cylToCart(I,r,x);

damageThreshold(I, probe_width, probe.lambda, "Propagated Probe:");
[f1,~] = computeDiameter(I, r, z, "Gauss","Propogated Probe Diameter");

% Checking for losses:
fprintf("\nUndisturbed Probe:");
[P_lost, T, P_back] = powerLoss(Iprobe, r);
fprintf("\nPropagated Probe:");
[Pp_lost, Tp, Pp_back] = powerLoss(I, r);

% Finding the z index where we get maximal intensity or focusing:
% [~, izMax,~] = findMax(Ipropogate_xz(:,:,itMax));
[~, izMax,~] = findMax(I_xz);
% prpFWHM = FWHM(I(:,:),r);
% [~,i] = min(prpFWHM);

PF_x_alongz(I_xz*1e-4,x,z,"Propagated Probe Intensity");
PF_x(I_xz(:,end)*1e-4,x,"Propagated Probe Intensity at the Sample End","Intensity [W/cm^2]");

%PF_x(Ipropogate_xz(:,izMax,itMax),x,"Probe Intensity at Maximum Focusing and 110[ps] delay");

%interactiveIxMovie(I_xz*1e-4, x, z, "Probe I(x,z)");
%PF_x(I_xz(:,i)*1e-4,x,"Probe Intensity at Maximum Focusing and 110[ps] delay","Intensity [W/cm^2]");
%PF_plot_xz(I_xz,z,x,"Probe Intensity 110[ps] delay",izMax);
% PF_x_alongz(Iprobe_xz,x,z,"Undisturbed Probe Intensity 110[ps] delay");

%nxz = cylToCart(n_complex(:,:,itMax),r,x);
%PF_x(imag(nxz(:,:)),x,"k 110[ps] delay");
%PF_x(real(nxz(:,:)),x,"n 110[ps] delay");

%PF_colormapAnimation_xz(Ipropogate*1e-4,r,z,x,1:Nt,"Probe Propogation", "Intensity [W/cm^2]", false, 50)

%PF_x(Ipropagate_xz(201,:,itMax),z,"Probe Intensity along x=0 at 110[ps] delay");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons and Other Tests
% =========================================================================

% Checking the best z0 and w0:
% z0_vec = (6:1:10)*1e-6;
% w0_vec = (1.17:0.005:1.2)*1e-6;
% CheckZ0(pump,z,r,x,t,probe,l);

% Checking the best pump energy:
% [opt_energy,opt_fwhm] = optimalEnergy(pump,z,r,x,t,probe);

% FWHM calculation:
%fwhm = zeros(Nz,Nt);

%for i = 1:Nt
%    fwhm(:,i) = FWHM(Ipropogate(:,:,i),r);
    %fwhm(:,i) = FWHM(Ipropagate_xz(:,:,i),x);
%end

%PF_FWHM_z(fwhm(:,itMax),z,sprintf("Propogated Probe FWHM at %d[ps]", t(itMax)*1e12));

% FWHM vs. t:
%tmin = min(fwhm(:));
% figure("Color",'w');
% plot(t,tmin(:))
% grid on; axis tight;
% xlabel("t [ps]"); ylabel("Minimal FWHM [\mum]");

% Comparing the probe with and without pump:
names = {"No Pump","35[nJ]"};
izList = {Nz,Nz};
PF_probeComparison({Iprobe_xz/max(Iprobe_xz(:,end)),I_xz/max(Iprobe_xz(:,end))}, names, x, z, t(itMax), izList);

figure('Color', 'w');
plot(z * 1e6, I(1, :) / max(Iprobe(1, :)), 'LineWidth', 2);
xlabel('z [\mum]');
ylabel('Intensity Ratio (I_{pump} / I_{no-pump})');
title('Probe Intensity Enhancement due to Pump');
grid on;

figure('Color', 'w');
plot(z * 1e6, I(1, :) ./ Iprobe(1, :), 'LineWidth', 2);
xlabel('z [\mum]');
ylabel('Intensity Ratio (I_{pump} / I_{no-pump})');
title('Probe Intensity Enhancement due to Pump');
grid on;

% Checking the maximal intensity vs. delay time:
%PF_maxIntensity(Ipropagate_xz, t);   


%PF_x(Ipropagate_xz(201,:,itMax),z,"Probe Intensity along x=0 at itmax")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pump Energy Related Effects
% =========================================================================

% Different pump energies:
% energies = [0, 10, 20, 30, 35, 64] * 1e-9;   % [J]
% names = {"0[nJ]","10[nJ]","20[nJ]","30[nJ]","35[nJ]","64[nJ]"};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

