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
dz = mean(diff(z));
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
pump_w0 = 1.329e-6;                     % m ~ 3; without BE ~ 3.345e-6
pump_z0 = 17e-6;%sp.Lz1                     % Beam waist location
l = 1;                                  % LG polynomial index

% Fresnel air-Si reflection:
% R = ((sp.n - 1)/(sp.n + 1))^2;
% pump_E = pump_E * (1 - R);

% Pump lasers:
pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, pump_w0, pump_z0, sp.n, l);  % 775[nm]

Ipump = pump.intensityProfileBLDumped(z);
Ipump_xz = cylToCart(Ipump,r,x);

% Pump diameters:
%computeDiameterTheoretically(z, "Donut",l, sp.wl2, pump_w0, z0, sp.n);
%computeDiameter(Ipump, r, z, "Donut");

% Checking we're under the destructive intensity:
damageThreshold(Ipump, pump_E, pump_width, pump.lambda, pump_w0, "Pump:", l);

% Plots:
%[~,f1] = PF_x(Ipump_xz(:,end)*1e-4,x,"Pump Intensity at Sample End","Intensity [W/cm^2]");
[~,f1] = PF_x(Ipump_xz(:,end)/max(Ipump_xz(:,end)),x,"Donut","","Intensity");
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
damageThreshold(Iund, probe_E, probe_width, probe.lambda, probe_w0, "Undisturbed Probe:");

% Plots:
[~,f2] = PF_x(Iund_xz(:,end)/max(Iund_xz(:,end)),x,"","Undisturbed Probe Intensity at the Sample End","Intensity [W/cm^2]");
%PF_x(Iund_xz(:,1)/max(Iund_xz(:,1)),x,"Undisturbed Probe Intensity at the Sample front","Intensity [W/cm^2]");
% PF_FWHM_z(FWHM(Iund,r),z,"Undisturbed Probe FWHM");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters Optimization and Other Tests
% =========================================================================

% Checking for the optimal waist, waist location and energy:
z0_vec = (0:1:20)*1e-6;%(-10:5:30)*1e-6;
w0_vec = (0.5:0.1:3)*1e-6;%(0.5:0.1:3.345)*1e-6;

%CheckZ0(pump,z,r,x,t,probe,z0_vec,l);
% CheckW0(pump,z,r,x,t,probe,w0_vec,l);
% CheckW0andZ0(pump,z,r,x,t,probe,z0_vec,w0_vec,0.91e-6,0.96e-6,"Figures_and_Results\Asaf\BE_sweep.txt");
%CheckW0andZ0(pump,z,r,x,t,probe,z0_vec,w0_vec,2.28e-6,1.27e-6,"Figures_and_Results\Asaf\pre_BE_sweep.txt");%0.91e-6,0.96e-6

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
% [~, ~, itMax] = findMax(pDiff);
t_delay = 100e-12;                      % Asaf's pump-probe delay time
[~, itMax] = min(abs(t - t_delay));
% pxz = cylToCart(pDiff(:,:,itMax),r,x) * 1e-6;

% [Nfcc, info] = intensityToCarriers(Ipump, pump_width, pump.lambda)ף
% Nfccxz = cylToCart(Nfcc,r,x) * 1e-6;

% PF_plot_xz(pxz,z,x,'pxz');
% PF_plot_xz(Nfccxz,z,x,'Nfcc');
% PF_plot_xz(Nfccxz./pxz,z,x);

n_complex = complexRefractiveIndex(pDiff(:,:,itMax), pump.lambda);%Nfcc

PF_complexRefractiveIndex(n_complex, r, z, x, pump.lambda, 1, t(itMax),true);

% PML mask parameters:
a = 4;
b = 8;

% Simulation for maximal delay (~90-110[ps]):
[~, Ipropagate] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex,a,b);
Ipropagate_xz = cylToCart(Ipropagate,r,x);

% Propagated probe Diameter:
%computeDiameter(Ipropagate, r, z, "Gauss");

% Checking we're under the destructive intensity:
damageThreshold(Ipropagate, probe_E, probe_width, probe.lambda, probe_w0, "Propagated Probe (l=1):");

% Checking for losses:
fprintf("\nUndisturbed Probe:");
[P_lost, T, P_back] = powerLoss(Iund, r);
fprintf("\nPropagated Probe:");
[Pp_lost, Tp, Pp_back] = powerLoss(Ipropagate, r);

[~, i] = findMax(Ipropagate_xz);        % Maximum intensity

% Plots:
[~,f3] = PF_x(Ipropagate_xz(:,end)/max(Iund_xz(:,end)),x,"Gauss","Propagated Probe Intensity at the Sample End","Intensity [W/cm^2]");
% PF_x_alongz(Ipropagate_xz*1e-4, x, z,"","Intensity [$\mathrm{W}/\mathrm{cm}^2$]");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons and Paper Plots
% =========================================================================

% Comparing different vortex orders:
% names = {"No Pump","l=1","l=2","l=3"};
% izList = {Nz,Nz,Nz,Nz};
% PF_probeComparison({Iund_xz/max(Iund_xz(:,end)),Ipropagate_xz/max(Iund_xz(:,end)), ...
%     I2(:,:,itMax)/max(Iund_xz(:,end)),I3(:,:,itMax)/max(Iund_xz(:,end))}, ...
%     names, x, z, t(itMax), izList);

% combined([hfig1, hfig2]);

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

function combined(figs)
    combinedFig = figure('Color','w');
    t = tiledlayout(combinedFig, 1, 2);
    
    titles = ["a)", "b)", "c)", "d)", "e)", "f)"];
    
    % Font controls
    tickFontSize  = 18;
    labelFontSize = 20;
    titleFontSize = 22;
    fontName = "Times New Roman";
    
    for k = 1:numel(figs)
    % Find the main axes in the source figure
    srcAx = findobj(figs(k), 'Type', 'Axes');
    srcAx = srcAx(1);

    % Create destination tile
    dstAx = nexttile(t, k);

    % Copy plotted content
    copyobj(allchild(srcAx), dstAx);

    srcLeg = findobj(figs(k), 'Type', 'Legend');
    if ~isempty(srcLeg)
        lgd = legend(dstAx, srcLeg(1).String, ...
            'Interpreter', 'latex', ...
            'Location', srcLeg(1).Location);
        lgd.FontSize = tickFontSize;
        lgd.FontName = fontName;
    end

    % Axis labels
    xlabel(dstAx, srcAx.XLabel.String, ...
        'Interpreter', 'latex', ...
        'FontSize', labelFontSize);

    ylabel(dstAx, srcAx.YLabel.String, ...
        'Interpreter', 'latex', ...
        'FontSize', labelFontSize);

    % Panel title: a), b), c), d)
    title(dstAx, titles(k), ...
        'Interpreter', 'latex', ...
        'FontSize', titleFontSize, ...
        'FontWeight', 'normal');

    % Copy important axis properties
    dstAx.XLim = srcAx.XLim;
    dstAx.YLim = srcAx.YLim;
    dstAx.XScale = srcAx.XScale;
    dstAx.YScale = srcAx.YScale;
    dstAx.YDir = srcAx.YDir;
    dstAx.Box = srcAx.Box;

    % Tick font controls
    dstAx.FontSize = tickFontSize;
    dstAx.FontName = fontName;
    dstAx.TickLabelInterpreter = 'latex';

    % Preserve aspect-ratio behavior
    dstAx.DataAspectRatio = srcAx.DataAspectRatio;
    dstAx.PlotBoxAspectRatio = srcAx.PlotBoxAspectRatio;
    dstAx.DataAspectRatioMode = srcAx.DataAspectRatioMode;
    dstAx.PlotBoxAspectRatioMode = srcAx.PlotBoxAspectRatioMode;

    % Copy colormap / CLim if relevant
    try
        dstAx.CLim = srcAx.CLim;
        colormap(dstAx, colormap(srcAx));
    end

    % Add colorbar if the original had one
    if ~isempty(findobj(figs(k), 'Type', 'ColorBar'))
        cb = colorbar(dstAx);
        cb.TickLabelInterpreter = 'latex';
        cb.FontSize = tickFontSize;
        cb.FontName = fontName;
    end

    % Make copied text annotations LaTeX too
    txt = findall(dstAx, 'Type', 'Text');
    for it = 1:numel(txt)
        txt(it).Interpreter = 'latex';
        txt(it).FontSize = tickFontSize;
        txt(it).FontName = fontName;
    end
end
end