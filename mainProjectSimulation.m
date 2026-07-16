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
tf = 1200;                              % 150 [ps]
t = (0:5:tf)*1e-12;                     % Time vector, [s]
Nt = numel(t);

% Spatial vectors intialization:
Nr = 1001;                              % Number of elements
Nz = 1001;

z = linspace(0, sp.Lz2, Nz);            % 25 micron

r = linspace(0,2*sp.Lx,Nr);             % 20x20 micron sample
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
pump_z0 = sp.Lz2;                       % Beam waist location
l = 1;                                  % LG polynomial index

% Fresnel air-Si reflection:
R = ((sp.n - 1)/(sp.n + 1))^2;
pump_E = pump_E * (1 - R);

% Pump laser beam:
pump = laser(sp.wl2, pump_width, pump_E, "Donut", r, phi, z, pump_w0, pump_z0, sp.n, l);  % 775[nm]

Ipump = pump.intensityProfileBLDumped(z);
Ipump_xz = cylToCart(Ipump,r,x);

% Plotting pump intensity:
%fig1 = PF_plot_xz(Ipump_xz/max(Ipump_xz(:)), z, x);
%[~,fig2] = PF_x(Ipump_xz(:,1)/max(Ipump_xz(:)), x,"Donut","Pump Intensity at the Sample Front","Intensity [au]");
%[~,fig3] = PF_x(Ipump_xz(:,end)/max(Ipump_xz(:,end)), x,"Donut","Pump Intensity at the Sample End","Intensity [au]");

%computeDiameter(Ipump, r, z, "Donut","Pump Diameters");
damageThreshold(Ipump, pump_E, pump_width, pump.lambda, pump_w0, "Pump:", l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FCC Creation & Diffusion
% =========================================================================

% Shining pump beam on the sample, causing e-h generation and diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)
% pDiff = intensityToCarriers(Ipump, pump_width, pump.lambda);    % Creates FCC distribution p(r,z,t)

% Finding the time when we get maximal concentration:
%[~, ~, itMax] = findMax(pDiff);         % tmax ~ 110[ps]
t_delay = 100e-12;                       % Nadav's pump-probe delay time
[~, itMax] = min(abs(t - t_delay));
%pxz = cylToCart(pDiff(:,:,itMax),r,x) * 1e-6;

%PF_x(pxz,x,"FCC Maximal Concentration in [1/cm^3] at t=110[ps]")
%PF_plot_xz(pxz,z,x,"FCC Maximal Concentration in [1/cm^3] at t=110[ps]");

%tIdx = [1, 19, 20, 50, 100, 150, 270, Nt];
%PF_colormapAnimation_xz(pDiff*1e-6, r, z ,x, 1:Nt, "FCC Concentration vs. Time", "FCC concentration [1/cm^3]", false, 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Refractive Index Calculation
% =========================================================================

% Calculating the complex refractive index n(r,z,t) (changes due to FCC generation):
n_complex = complexRefractiveIndex(pDiff(:,:,itMax), pump.lambda);

% Finding the z index where we get maximal absorption:
[~, izMax, ~] = findMax(imag(n_complex));%:,:,itMax)

%nxz = cylToCart(n_complex,r,x);
%PF_x(imag(nxz),x,"k 100[ps] delay");
%PF_x(real(nxz),x,"nr 100[ps] delay");

% PF_complexRefractiveIndex(n_complex, r, z, x, pump.lambda, izMax, t(itMax),true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probe Beam
% =========================================================================

% Probe parameters:
probe_width = 30e-12;                   % Pulse of 30[ps]
probe_E = pump.pulse_energy/100;        % [J]
probe_w0 = 2.267e-6;                    % [um], 2.378 
probe_z0 = sp.Lz2;

% Probe laser beam:
probe = laser(sp.wl2, probe_width, probe_E, "Gauss", r, phi, z, probe_w0, probe_z0, sp.n);   % 775[nm]

% Undisturbed probe propagating through the sample:
n = (sp.n + 1i*sp.alpha*probe.lambda/(4*pi))*ones(Nr,Nz);               % Includes BL 
[~, Iprobe] = propagationBPM_rz(probe.profile(:,1),r,z,probe,n,0,8);    % [W/m^2]
Iprobe_xz = cylToCart(Iprobe,r,x);

% L = (5:5:135)*1e-9;
% 
% a = 4;
% b = 8;
% t_delay = 100e-12;
% [~, itDelay] = min(abs(t - t_delay));
% fwhm=zeros(1,numel(L)+1);
% 
% colors = [sp.colors(12), sp.colors(13), sp.colors(17)];
% for i = 1:numel(L)
% 
%     fprintf(' %.3f\n', L(i)*1e6);
% 
%     pump = laser(sp.wl2, pump_width, L(i), ...
%         "Donut", r, phi, z, pump_w0, pump_z0, sp.n, l);%z0_best(i)
%     Ipump = pump.intensityProfileBLDumped(z);
% 
% 
%     probe = laser(sp.wl2, probe_width, L(i)/100, ...
%         "Gauss", r, phi, z, probe_w0, probe_z0, sp.n);
%     % Carrier distribution
%     pDiff = FCCDiffusion(pump, t, r, z);%intensityToCarriers(Ipump, pump_width, pump.lambda);
%     n_complex = complexRefractiveIndex(pDiff(:,:,itDelay), pump.lambda);
% 
%     % [Nfcc, ~] = FCC_fromGaussianPulse(Ipump);
%     % n_complex = complexRefractiveIndex(Nfcc, pump.lambda);%pDiff(:,:,itMax)
% 
%     % Probe propagation after pump
%     [~, I_afterPump] = propagationBPM_rz( ...
%         probe.profile(:,1), r, z, probe, n_complex, a, b);
%     I_after_xz = cylToCart(I_afterPump,r,x);
% 
%     %temp = FWHM(I_afterPump,r);
%     % fwhm(i)=temp(end);
%     [fwhm(i+1),~]=PF_x(I_after_xz(:,end)*1e4,x,"Propagated Probe Intensity at the Sample End","Intensity [W/cm^2]");
% 
% end
% hfig2 = figure('Color','w');
% 
% plot([0,L*1e9], fwhm*1e6,'LineWidth', 4,'Color', colors(2));
% 
% grid on; axis tight;
% 
% xlabel("Pump Energy [nJ]", 'Interpreter','latex');
% ylabel("Probe FWHM $[\mu\mathrm{m}]$", 'Interpreter','latex');
% 
% picturewidth = 25;
% hw_ratio = 0.65;
% 
% set(findall(hfig2,'-property','FontSize'), 'FontSize', 20);
% set(findall(hfig2,'-property','TickLabelInterpreter'), ...
%     'TickLabelInterpreter','latex');
% set(hfig2, 'Units','centimeters', ...
%     'Position',[3 3 picturewidth hw_ratio*picturewidth]);

fig4 = PF_plot_xz(Iprobe_xz/max(Iprobe_xz(:)), z, x);
% [~,fig5] = PF_x(Iprobe_xz(:,1)/max(Iprobe_xz(:)), x,'',"Undisturbed Probe Intensity at the Sample Front","Intensity [au]");
[~,fig6] = PF_x(Iprobe_xz(:,end)/max(Iprobe_xz(:,end)), x,'',"Undisturbed Probe Intensity at the Sample End","Intensity [au]");
%PF_x_alongz(Iprobe_xz*1e-4,x,z,"Propagated Probe Intensity");

%absorptionDepth(Iprobe(1,:),z,probe.lambda);
damageThreshold(Iprobe, probe_E, probe_width, probe.lambda, probe_w0, "Undisturbed Probe:");

% Probe Diameter:
%computeDiameterTheoretically(z, "Gauss",0, sp.wl2, probe_w0, 0, sp.n);
%computeDiameter(Iprobe, r, z, "Gauss");
%PF_FWHM_z(FWHM(Iprobe,r),z,"Undisturbed Probe FWHM");

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
[~, I] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex,a,b);%(:,:,itMax)
I_xz = cylToCart(I,r,x);

damageThreshold(I, probe_E, probe_width, probe.lambda, probe_w0, "Propagated Probe:");
[fig7,~] = computeDiameter(I, r, z, "Gauss","Propogated Probe Diameter");

fig8 = PF_plot_xz(I_xz/max(Iprobe_xz(:)), z, x); 
[~,fig9] = PF_x(I_xz(:,end)/max(I_xz(:,end)), x,'',"Propagated Probe Intensity at the Sample End","Intensity [au]");

% Checking for losses:
fprintf("\nUndisturbed Probe:");
[P_lost, T, P_back] = powerLoss(Iprobe, r);
fprintf("\nPropagated Probe:");
[Pp_lost, Tp, Pp_back] = powerLoss(I, r);

% Finding the z index where we get maximal intensity or focusing:
% [~, izMax,~] = findMax(Ipropogate_xz(:,:,itMax));
[~, izMax,~] = findMax(I_xz);

PF_x_alongz(I_xz*1e-4,x,z,"Propagated Probe Intensity");

%PF_x(Ipropogate_xz(:,izMax,itMax),x,"Probe Intensity at Maximum Focusing and 110[ps] delay");

%interactiveIxMovie(I_xz*1e-4, x, z, "Probe I(x,z)");
%PF_x(I_xz(:,i)*1e-4,x,"Probe Intensity at Maximum Focusing and 110[ps] delay","Intensity [W/cm^2]");
%PF_plot_xz(I_xz,z,x,"Probe Intensity 110[ps] delay",izMax);
% PF_x_alongz(Iprobe_xz,x,z,"Undisturbed Probe Intensity 110[ps] delay");

%PF_colormapAnimation_xz(I*1e-4,r,z,x,1:Nt,"Probe Propogation", "Intensity [W/cm^2]", false, 25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparisons and Other Tests
% =========================================================================

% Checking the best z0 and w0:
% z0_vec = (0:1:25)*1e-6;
% w0_vec = (1.17:0.005:1.2)*1e-6;
% CheckZ0(pump,z,r,x,t,probe,l);

% Checking the best pump energy:
% [opt_energy,opt_fwhm] = optimalEnergy(pump,z,r,x,t,probe);

% Comparing the probe with and without pump:
% names = {"No Pump","35[nJ]"};
% izList = {Nz,Nz};
% fig10 = PF_probeComparison({Iprobe_xz/max(Iprobe_xz(:,end)),I_xz/max(Iprobe_xz(:,end))}, names, x, z, t(itMax), izList,"Intensity [au]");

% combinedFig = figure('Color','w');
% t = tiledlayout(combinedFig, 2, 2);
% 
% figs = [fig8, fig9, fig7, fig10];
% titles = ["a)", "b)", "c)", "d)", "e)", "f)"];
% 
% % Font controls
% tickFontSize  = 18;
% labelFontSize = 20;
% titleFontSize = 22;
% fontName      = "Times New Roman";
% 
% for k = 1:numel(figs)
%     % Find the main axes in the source figure
%     srcAx = findobj(figs(k), 'Type', 'Axes');
%     srcAx = srcAx(1);
% 
%     % Create destination tile
%     dstAx = nexttile(t, k);
% 
%     % Copy plotted content
%     copyobj(allchild(srcAx), dstAx);
% 
%     srcLeg = findobj(figs(k), 'Type', 'Legend');
%     if ~isempty(srcLeg)
%         lgd = legend(dstAx, srcLeg(1).String, ...
%             'Interpreter', 'latex', ...
%             'Location', srcLeg(1).Location);
%         lgd.FontSize = tickFontSize;
%         lgd.FontName = fontName;
%     end
% 
%     % Axis labels
%     xlabel(dstAx, srcAx.XLabel.String, ...
%         'Interpreter', 'latex', ...
%         'FontSize', labelFontSize);
% 
%     ylabel(dstAx, srcAx.YLabel.String, ...
%         'Interpreter', 'latex', ...
%         'FontSize', labelFontSize);
% 
%     % Panel title: a), b), c), d)
%     title(dstAx, titles(k), ...
%         'Interpreter', 'latex', ...
%         'FontSize', titleFontSize, ...
%         'FontWeight', 'normal');
% 
%     % Copy important axis properties
%     dstAx.XLim = srcAx.XLim;
%     dstAx.YLim = srcAx.YLim;
%     dstAx.XScale = srcAx.XScale;
%     dstAx.YScale = srcAx.YScale;
%     dstAx.YDir = srcAx.YDir;
%     dstAx.Box = srcAx.Box;
% 
%     % Tick font controls
%     dstAx.FontSize = tickFontSize;
%     dstAx.FontName = fontName;
%     dstAx.TickLabelInterpreter = 'latex';
% 
%     % Preserve aspect-ratio behavior
%     dstAx.DataAspectRatio = srcAx.DataAspectRatio;
%     dstAx.PlotBoxAspectRatio = srcAx.PlotBoxAspectRatio;
%     dstAx.DataAspectRatioMode = srcAx.DataAspectRatioMode;
%     dstAx.PlotBoxAspectRatioMode = srcAx.PlotBoxAspectRatioMode;
% 
%     % Copy colormap / CLim if relevant
%     try
%         dstAx.CLim = srcAx.CLim;
%         colormap(dstAx, colormap(srcAx));
%     end
% 
%     % Add colorbar if the original had one
%     if ~isempty(findobj(figs(k), 'Type', 'ColorBar'))
%         cb = colorbar(dstAx);
%         cb.TickLabelInterpreter = 'latex';
%         cb.FontSize = tickFontSize;
%         cb.FontName = fontName;
%     end
% 
%     % Make copied text annotations LaTeX too
%     txt = findall(dstAx, 'Type', 'Text');
%     for it = 1:numel(txt)
%         txt(it).Interpreter = 'latex';
%         txt(it).FontSize = tickFontSize;
%         txt(it).FontName = fontName;
%     end
% end
% 

% figure('Color', 'w');
% plot(z * 1e6, I(1, :) / max(Iprobe(1, :)), 'LineWidth', 2);
% xlabel('z [\mum]');
% ylabel('Intensity Ratio (I_{pump} / I_{no-pump})');
% title('Probe Intensity Enhancement due to Pump');
% grid on;
% 
% figure('Color', 'w');
% plot(z * 1e6, I(1, :) ./ Iprobe(1, :), 'LineWidth', 2);
% xlabel('z [\mum]');
% ylabel('Intensity Ratio (I_{pump} / I_{no-pump})');
% title('Probe Intensity Enhancement due to Pump');
% grid on;

