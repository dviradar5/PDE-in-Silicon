%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation for Asaf's Paper
% -------------------------------------------------------------------------
% In this file we run the simulation with the parameters from Asaf's paper,
% with the goal of seeing good correlation between the results
% 
% We assume the Beam Expender was deployed and use the parameter values
% that it created. In order to implement other BE orders into the
% simulation, we simply devide the original pump w0 by the order we want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters and Vector Definitions
% =========================================================================

sp = systemParameters();

% Temporal vector intialization:
tf = 200;                               % 200[ps]
t = (0:5:tf)*1e-12;                     % Time vector, [s]
Nt = numel(t);

% Spatial vectors intialization:
Nr = 1001;                              % Number of elements
Nz = 1001;

z = linspace(0, sp.Lz1, Nz);            % 20 micron

r = linspace(0,sp.Lx,Nr);               % 10x10 micron sample
phi = atan(1);                          % y = x, radial symmetry

x = linspace(-sp.Lx, sp.Lx, 2*Nr-1);
Nx = numel(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pump Beam
% =========================================================================

% Pump parameters:
pump_width = 30e-12;                    % Pulse of 30[ps]
pump_E = 20e-9;                         % Beam total energy of 20[nJ], [J]
pump_w0 = 3.345e-6;                     % m ~ 3; without BE ~ 1.329e-6   
pump_z0 = 17e-6;                        % Beam waist location
l = 1;                                  % LG polynomial index

% Fresnel air-Si reflection:
%R = ((sp.n - 1)/(sp.n + 1))^2;
%pump_E = pump_E * (1 - R);

% Probe parameters:
probe_width = 30e-12;                   % Pulse of 50[ps]
probe_w0 = 5e-6;                        % 1.875e-6, 10e-6          
probe_z0 = sp.Lz1;
n = (sp.n + 1i*sp.alpha*sp.wl2/(4*pi))*ones(Nr,Nz);     % Includes BL 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ckptFile = 'ckpt_asaf_probe_FWHM_vs_pump_diameter.mat';
resetCheckpoint = false;   % set true if you want to start from scratch

if resetCheckpoint && isfile(ckptFile)
    delete(ckptFile);
    fprintf('Deleted old checkpoint file: %s\n', ckptFile);
end

% Different energies:
%L = [0,50e-9,75e-9,100e-9];
%labels = ["No Pump","50 [nJ]","75 [nJ]","100 [nJ]"];

% Different pump radii:
L = (0.5:0.25:8)*1e-6;%[0,1.185e-6,3.385e-6,5e-6,7.5e-6];%,1.185e-6
%labels = ["No Pump","1.185 [\mum]","3.385 [\mum]","5 [\mum]","7.5 [\mum]"];%,"1.185 [\mum]"

% Different topological charge:
%L = [1,2];
%labels = ['1','2'];

Dfigs = gobjects(1,numel(L));
Ifigs = gobjects(1,numel(L));

w0_vec = [1.875e-6,5e-6,10e-6];
%end_fwhm = zeros(numel(w0_vec),numel(L));

if isfile(ckptFile)
    load(ckptFile, ...
        'L', 'w0_vec', 'end_fwhm', 'completed', ...
        'D', 'l', 'pump_E', 'pump_z0', 'probe_z0');

    fprintf('Loaded checkpoint file: %s\n', ckptFile);

    % Safety checks
    if ~isequal(L, (0.5:0.25:8)*1e-6)
        warning('Loaded checkpoint uses a different L vector than the current script.');
    end

    if ~isequal(w0_vec, [1.875e-6, 5e-6, 10e-6])
        warning('Loaded checkpoint uses a different w0_vec than the current script.');
    end

else
    end_fwhm = NaN(numel(w0_vec), numel(L));
    completed = false(numel(w0_vec), numel(L));
    D = NaN(size(L));  % filled later

    save(ckptFile, ...
        'L', 'w0_vec', 'end_fwhm', 'completed', ...
        'D', 'l', 'pump_E', 'pump_z0', 'probe_z0');

    fprintf('Created new checkpoint file: %s\n', ckptFile);
end

for j = 1:numel(w0_vec)
    probe = laser(sp.wl2, probe_width, pump_E/100, "Gauss", r, phi, z, w0_vec(j), probe_z0, sp.n);   % 775[nm]

    for i = 1:numel(L)
        %pump = laser(sp.wl2, pump_width, L(i), "Donut", r, phi, z, pump_w0, z0, sp.n, l);
        %pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, pump_w0, z0, sp.n, L(i));

        pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, L(i), pump_z0, sp.n, l);

        %Ipump = pump.intensityProfileBLDumped(z);
        %Ipump_xz = cylToCart(Ipump,r,x);
        %PF_plot_xz(Ipump_xz/max(Ipump_xz(:)), z, x, "Pump Normalized");

        pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)
        [~, ~, itMax] = findMax(pDiff); 

        n_complex = complexRefractiveIndex(pDiff, pump.lambda);
        %PF_complexRefractiveIndex(n_complex(:,:,itMax), r, z, x, pump.lambda, 1, t(itMax));

        if L(i) == 0
            %probe = laser(sp.wl2, probe_width, 1e-9, "Gauss", r, phi, z, probe_w0, 0, sp.n);
            [~, I] = propagationBPM_rz(probe.profile(:,1),r,z,probe,n,0,8);    % [W/m^2]
        else
            %probe = laser(sp.wl2, probe_width, L(i)/100, "Gauss", r, phi, z, probe_w0, 0, sp.n);
            [~, I] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,itMax),4,20);
        end

        %PF_plot_xz(I_xz*1e-4, z, x, "Probe Intensity");

        %Imax = max(I_xz);% .* exp(sp.alpha * z);
        % Ifigs(i) = figure('Color','w');
        % plot(z*1e6,Imax/Imax(1),'LineWidth', 3);
        % xlabel("z [\mum]"); ylabel("Intensity [au]");
        % title(labels(i));

        temp = FWHM(I,r, 5*mean(diff(r)), 0.02);
        end_fwhm(j,i) = temp(end);
        %PF_x_alongz(I_xz*1e-4,x,z,"Propagated Probe Intensity");

        %[Dfigs(i),~] = computeDiameter(I, r, z, "Gauss","Propogated Probe Diameter");

        completed(j,i) = true;

        % Radius to LG diameter:
        D = 2*sqrt(abs(l)/2) .* L;

        save(ckptFile, ...
            'L', 'w0_vec', 'end_fwhm', 'completed', ...
            'D', 'l', 'pump_E', 'pump_z0', 'probe_z0');

        fprintf('Saved checkpoint after j = %d, i = %d\n', j, i);
    end
end
PF_combinePlots(L, labels, Dfigs,"Width [\mum]" , '');
PF_combinePlots(L, labels, Ifigs,"Intensity [au]" , '');

ckptFile = 'ckpt_asaf_probe_FWHM_vs_pump_diameter.mat';
load(ckptFile, 'D', 'end_fwhm', 'w0_vec');

hfig = figure('Color','w');

plot(D*1e6, end_fwhm(1,:)*1e6, 'r-',  'LineWidth', 4); hold on;
plot(D*1e6, end_fwhm(2,:)*1e6, 'g--', 'LineWidth', 4);
plot(D*1e6, end_fwhm(3,:)*1e6, 'b:',  'LineWidth', 4); hold off;

xlabel("pump Diameter [$\mu$m]", 'Interpreter','latex');
ylabel("probe FWHM [$\mu$m]", 'Interpreter','latex');

legend( ...
    '$\omega_{0,\mathrm{probe}} = 1.875\,[\mu\mathrm{m}]$', ...
    '$\omega_{0,\mathrm{probe}} = 5\,[\mu\mathrm{m}]$', ...
    '$\omega_{0,\mathrm{probe}} = 10\,[\mu\mathrm{m}]$', ...
    'Interpreter','latex', ...
    'Location','best');

grid on; axis tight; %box off;

fname = 'asaf_probe_FWHM_vs_pump_diameter';

picturewidth = 25;
hw_ratio = 0.65;

set(findall(hfig,'-property','FontSize'), 'FontSize', 20);
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(hfig, 'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);

% Better modern export
exportgraphics(hfig, "Figures_and_Results\" + fname + ".png", 'Resolution', 300);

% Pump lasers:
pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, pump_w0, pump_z0, sp.n, l);  % 775[nm]

Ipump = pump.intensityProfileBLDumped(z);
Ipump_xz = cylToCart(Ipump,r,x);

% Pump diameters:
%computeDiameterTheoretically(z, "Donut",l, sp.wl2, pump_w0, z0, sp.n);
%computeDiameter(Ipump, r, z, "Donut");

% Checking we're under the destructive intensity:
damageThreshold(Ipump, pump_width, pump.lambda, "Pump:");

% Plots:
%[~,f1] = PF_x(Ipump_xz(:,end)*1e-4,x,"Pump Intensity at Sample End","Intensity [W/cm^2]");
[~,f1] = PF_x(Ipump_xz(:,end)/max(Ipump_xz(:,end)),x,"","Intensity");
%PF_x_alongz(Ipump_xz*1e-4, x, z,"Pump Intensity","Intensity [W/cm^2]");
%PF_FWHM_z(FWHM(Ipump,r),z)                                         %FIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probe Beam
% =========================================================================

% Probe parameters:
probe_width = 30e-12;                   % Pulse of 30[ps]
probe_E = pump.pulse_energy/100;        % [J]
probe_w0 = 1.875e-6;      
probe_z0 = sp.Lz1;

% Probe laser:
probe = laser(sp.wl2, probe_width, probe_E, "Gauss", r, phi, z, probe_w0, probe_z0, sp.n);   % 775[nm]

% Undisturbed probe propagating through the sample:
n = (sp.n + 1i*sp.alpha*sp.wl2/(4*pi))*ones(Nr,Nz);                 % Includes BL 
[~, Iund] = propagationBPM_rz(probe.profile(:,1),r,z,probe,n,4,8);  % [W/m^2]
Iund_xz = cylToCart(Iund,r,x);

% Checking we're under the destructive intensity:
damageThreshold(Iund, probe_width, probe.lambda, "Undisturbed Probe:");

% Plots:
[~,f2] = PF_x(Iund_xz(:,end)/max(Iund_xz(:,end)),x,"Undisturbed Probe Intensity at the Sample End","Intensity [W/cm^2]");
%PF_x(Iund_xz(:,1)/max(Iund_xz(:,1)),x,"Undisturbed Probe Intensity at the Sample front","Intensity [W/cm^2]");
% PF_FWHM_z(FWHM(Iund,r),z,"Undisturbed Probe FWHM");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters Optimization and Other Tests
% =========================================================================

% Checking for the optimal waist, waist location and energy:
z0_vec = (0:1:20)*1e-6;%(-10:5:30)*1e-6;
w0_vec = (3:0.01:4)*1e-6;%(0.5:0.1:3.345)*1e-6;

%CheckZ0(pump,z,r,x,t,probe,z0_vec,l);
% CheckW0(pump,z,r,x,t,probe,w0_vec,l);
%CheckW0andZ0(pump,z,r,x,t,probe,z0_vec,w0_vec,0.91e-6,0.96e-6,"Figures_and_Results\BE_sweep.txt");
%CheckW0andZ0(pump,z,r,x,t,probe,z0_vec,w0_vec,2.28e-6,1.27e-6,"Figures_and_Results\pre_BE_sweep.txt");%0.91e-6,0.96e-6

% [opt_energy,opt_fwhm] = optimalEnergy(pump,z,r,x,t,probe);

% Checking different pump vortex orders:
%l_vec = 1:1:8;
%CheckLGOrder(pump,z,r,x,t,probe,l_vec);

% Probe Diameter:
%computeDiameterTheoretically(z, "Gauss",0, sp.wl2, probe_w0, 0, sp.n);
%computeDiameter(Iund, r, z, "Gauss");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Run
% =========================================================================

% FCC diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)
[~, ~, itMax] = findMax(pDiff);

t_delay = 100e-12;                      % Asaf's pump-probe delay time
[~, itMax] = min(abs(t - t_delay));

% Calculating the new complex refractive index:
n_complex = complexRefractiveIndex(pDiff, pump.lambda);
[~, in] = findMax(n_complex(:,:,itMax));% Maximum intensity

% PML mask parameters:
a = 4;
b = 8;

% Simulation for maximal delay (~90-110[ps]):
[~, Ipropagate] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,itMax),a,b);
Ipropagate_xz = cylToCart(Ipropagate,r,x);

% Propagated probe Diameter:
%computeDiameter(Ipropagate, r, z, "Gauss");

% Checking we're under the destructive intensity:
damageThreshold(Ipropagate, probe_width, probe.lambda, "Propagated Probe (l=1):");

% Checking for losses:
fprintf("\nUndisturbed Probe:");
[P_lost, T, P_back] = powerLoss(Iund, r);
fprintf("\nPropagated Probe:");
[Pp_lost, Tp, Pp_back] = powerLoss(Ipropagate, r);

[~, i] = findMax(Ipropagate_xz);        % Maximum intensity

% Plots:
[~,f3] = PF_x(Ipropagate_xz(:,end)/max(Iund_xz(:,end)),x,"Propagated Probe Intensity at the Sample End","Intensity [W/cm^2]");
PF_x_alongz(Ipropagate_xz*1e-4, x, z,"Propagated Probe Intensity","Intensity [W/cm^2]");

% prpFWHM = FWHM(Ipropagate,r);
%[~,i] = min(prpFWHM);                  % Minimal FWHM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons and Paper Plots
% =========================================================================

% Comparing different vortex orders:
% names = {"No Pump","l=1","l=2","l=3"};
% izList = {Nz,Nz,Nz,Nz};
% PF_probeComparison({Iund_xz/max(Iund_xz(:,end)),Ipropagate_xz/max(Iund_xz(:,end)), ...
%     I2(:,:,itMax)/max(Iund_xz(:,end)),I3(:,:,itMax)/max(Iund_xz(:,end))}, ...
%     names, x, z, t(itMax), izList);

% Combining the figures:                    MAKE PF #######################
% figure("Color",'w');
% ax1 = subplot(1,3,1);
% copyobj(allchild(findobj(f2,"Type","axes")), ax1);
% ylabel("Normalized Intensity", "FontSize",15); ylim([0,1.1]);
% xlabel("x [\mum]", "FontSize",15);
% 
% ax2 = subplot(1,3,2);
% copyobj(allchild(findobj(f1,"Type","axes")), ax2);
% ylabel("Normalized Intensity", "FontSize",15); ylim([0,1.1]);
% xlabel("x [\mum]", "FontSize",15);
% 
% ax3 = subplot(1,3,3);
% copyobj(allchild(findobj(f3,"Type","axes")), ax3);
% ylabel("Normalized Intensity", "FontSize",15); ylim([0,1.1]);
% xlabel("x [\mum]", "FontSize",15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other Plots
% =========================================================================

% Pump:
PF_x_alongz(Ipump_xz*1e-4, x, z,"Pump Intensity","Intensity [W/cm^2]");
PF_x(Ipump_xz(:,1)*1e-4,x,"Pump Intensity at the Sample Surface (z=0)", "Intensity [W/cm^2]");
PF_beamIntensityPlot_xy(Ipump(:,1)*1e-4, r, x);

% Original probe:
PF_x_alongz(Iund_xz*1e-4, x, z,"Undisturbed Probe Intensity","Intensity [W/cm^2]");
PF_x(Iund_xz(:,1)*1e-4,x,"Undisturbed Probe Intensity at the Sample Surface (z=0)","Intensity [W/cm^2]");
PF_x(Iund_xz(:,i)*1e-4,x,"Undisturbed Probe Intensity at Maximum Intensity","Intensity [W/cm^2]",[]);
PF_beamIntensityPlot_xy(Iund(:,1)*1e-4, r, x,"Undisturbed Probe Intensity at the Sample Surface (z=0) [W/cm^2]");

% Complex refractive index:
PF_complexRefractiveIndex(n_complex(:,:,itMax), r, z, x, pump.lambda, 1, t(itMax));

% Propagated probe:
PF_x_alongz(Ipropagate_xz*1e-4, x, z,"Propagated Probe Intensity","Intensity [W/cm^2]");
PF_x(Ipropagate_xz(:,i)*1e-4,x,"Propagated Probe Intensity at Maximum Intensity","Intensity [W/cm^2]");
PF_x(Ipropagate_xz(:,Nz)*1e-4,x,"Propagated Probe Intensity at the Sample End","Intensity [W/cm^2]");
PF_beamIntensityPlot_xy(Ipropagate(:,i)*1e-4, r, x, "Propagated Probe Intensity at Maximum Intensity [W/cm^2]");

% FWHM:
PF_FWHM_z(prpFWHM,z,"Propagated Probe FWHM");
