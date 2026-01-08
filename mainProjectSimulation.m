clc; clear;

sp = systemParameters();

% Temporal vector intialization:
tf = 30000;                 % 30 [ns]
t = 0:10:tf;                % time vector [ps]
t = t .* 1e-12;             % time vector [s]
Nt = numel(t);

% Spacial vector intialization:
Nr = 201;                               % Number of steps
Nz = 201;

z = linspace(0, sp.Lz, Nz);
r = linspace(0,sp.Lx/sqrt(2),Nr);       % FIX
phi = atan(1);                          % y = x due to radial symmetry

% Laser constants:
pmp_wd = 3e-11;                         % pulse of 30[ps]
pmp_D = sp.D2;
w0 = 5e-6;
pump_E0 = sqrt(sp.pulse_energies(5));   % E0 = 1

% Pump lasers:
%pump1 = laser(sp.wl1, pmp_wd, pmp_D, "Donut", r, phi, z, w0, 0, pump_E0);  % z0 = 0
pump2 = laser(sp.wl2, pmp_wd, pmp_D, "Donut", r, phi, z, w0, 0, pump_E0);  % 775[nm]
%pump3 = laser(sp.wl2, pmp_wd, pmp_D, "Gauss", r, phi, z, w0, 0, pump_E0);

% Shining 775[nm] donut beam on the sample:
p2diff = FCCDiffusion(pump2, t, r, z);

% iz = 1;         % z(1) = 0
% it = 1;         % t(1) = 0
% fprintf('z(iz) = %d', z(iz));

% FCC distribution at different times:
% figure;
% subplot(1,3,1);
% [P, xPlot, yPlot] = cylToCart(p2diff(:,:,1), r, iz);
% imagesc(xPlot, yPlot, P);
% colorbar;
% subplot(1,3,2);
% [P, xPlot, yPlot] = cylToCart(p2diff(:,:,(Nt+1)/2), r, iz);
% imagesc(xPlot, yPlot, P);
% colorbar;
% subplot(1,3,3);
% [P, xPlot, yPlot] = cylToCart(p2diff(:,:,Nt), r, iz);
% imagesc(xPlot, yPlot, P);
% colorbar;

% Comparing different times and with and without Beer-Lambert dumping:
% P1 = cylToCart(p2diff(:,:,1), r, iz0, xPlot, yPlot);
% P2 = cylToCart(p2diff(:,:,Nt), r, iz0, xPlot, yPlot);
% 
% I1 = abs(P1).^2;
% I2 = abs(P2).^2;
% 
% % intensity dumping according to Beer-Lambert's law:
% I1BL = I1 * exp(-z(iz0)/sp.penDepth);
% I2BL = I2 * exp(-z(iz0)/sp.penDepth);
% 
% figure;
% subplot(2,1,1);
% imagesc(xPlot*1e6, yPlot*1e6, I1);
% axis image; set(gca,'YDir','normal');
% xlabel('x [\mum]'); ylabel('y [\mum]');
% title('p(x,y) at z=0');
% colorbar;
% 
% subplot(2,1,2);
% imagesc(xPlot*1e6, yPlot*1e6, I1BL);
% axis image; set(gca,'YDir','normal');
% xlabel('x [\mum]'); ylabel('y [\mum]');
% title('p(x,y) at z=0 dumped');
% colorbar;

% Probe laser:
prb_wd = 5e-11;   % pulse of 50[ps]
prb_D = sp.D1;
probe_E0 = sqrt(sp.pulse_energies(5));  % E0=1
probe = laser(sp.wl1, prb_wd, prb_D, "Gauss", r, phi, z, w0, 0, probe_E0);  % 1550nm

% Shining the probe on the sample:
prb_prf_after = probeAbsorptionBS(probe, p2diff, r, z, t);

[P1, xPlot1, yPlot1] = cylToCart(prb_prf_after(:,:,1), r, 1);
I1 = (abs(P1).^2) .* exp(-z(1)/sp.penDepth);
[P2, xPlot2, yPlot2] = cylToCart(prb_prf_after(:,:,(Nt+1)/2), r, 1);
I2 = (abs(P2).^2) .* exp(-z(1)/sp.penDepth);
figure;
subplot(1,2,1);
surf(xPlot1*1e6, yPlot1*1e6, I1);
shading interp
xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [au]')
colorbar
subplot(1,2,2);
surf(xPlot2*1e6, yPlot2*1e6, I2);
shading interp
xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [au]')
colorbar

% iz=1;
% 
% [P1, xPlot1, yPlot1] = cylToCart(prb_prf_after(:,:,1), r, iz);
% [P2, xPlot2, yPlot2] = cylToCart(prb_prf_after(:,:,Nt), r, iz);
% I1 = P1 * exp(-z(iz)/sp.penDepth);
% I2 = P2 * exp(-z(iz)/sp.penDepth);
% 
% figure;
% subplot(1,2,1);
% surf(xPlot1*1e6, yPlot1*1e6, I1);
% shading interp
% xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [J]')
% title('probe at t=0');
% colorbar
% 
% subplot(1,2,2);
% surf(xPlot2*1e6, yPlot2*1e6, I2);
% shading interp
% xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [J]')
% title('probe at t=30[ns]');
% colorbar

% % Comparing with and without pump:
% [P1, xPlot1, yPlot1] = cylToCart(probe.profile, r, iz);
% [P2, xPlot2, yPlot2] = cylToCart(prb_prf_after(:,:,it), r, iz);
% I1 = abs(P1).^2;
% I2 = abs(P2).^2;
% I1BL = I1 * exp(-z(iz)/sp.penDepth);
% I2BL = I2 * exp(-z(iz)/sp.penDepth);
% 
% figure;
% subplot(1,2,1);
% surf(xPlot1*1e6, yPlot1*1e6, I1BL);
% shading interp
% xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [J]')
% title('probe at z=0 without pump');
% colorbar
% 
% subplot(1,2,2);
% surf(xPlot2*1e6, yPlot2*1e6, I2BL);
% shading interp
% xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [J]')
% title('probe at z=0 with pump');
% colorbar

% FWHM:
% I1_cyl = abs(probe.profile).^2;
% I2_cyl = abs(prb_prf_after(:,:,it)).^2;
% I1cylBL = I1_cyl * exp(-z(iz)/sp.penDepth);
% I2cylBL = I2_cyl * exp(-z(iz)/sp.penDepth);
% 
% fwhm1 = FWHM(I1cylBL, r);
% fwhm2 = FWHM(I2cylBL, r);
% 
% figure;
% plot(z*1e3, fwhm1*1e6, 'g');
% hold on;
% plot(z*1e3, fwhm2*1e6, 'm');
% hold off;
% xlabel('z [mm]'); ylabel('FWHM [\mum]');
% legend('without pump','with pump');
% title('FWHM of probe with and without pump');
% grid on;

% pump-probe delay:
delay = 0:5:25;
delay = delay * 1e-9;   % delay in [ns]

zi = 41; % 5 micron into the sample

% figure;
% for i = 1:length(delay)
%     [P, xPlot, yPlot] = cylToCart(prb_prf_after(:,:,(300*i+1)), r, zi);
%     I = (abs(P).^2) .* exp(-z(zi)/sp.penDepth);
% 
%     subplot(2,length(delay)/2,i);
%     surf(xPlot*1e6, yPlot*1e6, I);
%     shading interp
%     xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [au]')
%     title(sprintf('\\DeltaT = %.2d[s]',delay(i)));
%     colorbar
% end
% 
% sgtitle(sprintf('Pump-Probe delay effect on probe intensity at z=%.2d[\\mum]', z(zi)*1e6));
