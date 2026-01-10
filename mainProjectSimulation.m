clc; clear;

sp = systemParameters();

% Temporal vector intialization:
tf = 30000;                 % 0-30 [ns]
t = 0:10:tf;                % time vector [ps]
t = t .* 1e-12;             % time vector [s]
dt = t(2)-t(1);
Nt = numel(t);

% Spacial vector intialization:
Nr = 201;                               % Number of elements
Nz = 201;

z = linspace(0, sp.Lz, Nz);
r = linspace(0,sqrt(2)*sp.Lx,Nr);       % Sample is 10x10 micron, so this r covers it
phi = atan(1);                          % y = x, radial symmetry
x = linspace(-max(r), max(r), 2*Nr-1);

% Source:
y_src = -1;
polarization = Axis.x;
src =  RectSrc(Axis.y, y_src, [0 0.02; -0.1 0.1], polarization);

% Laser constants:
pmp_wd = 3e-11;        % pulse of 30[ps]
pmp_D = sp.D2;
w0 = 5e-6;
z0 = -7;        % Beam waist location
pump_E = 16e-9;   % [J] FIXXXXXXXXXXXXX
dE = 8e-9;      % Size of each energy portion [J]
nsim = round(pump_E/dE);    % Number of iterations

% Pump laser:
pump = laser(sp.wl2, pmp_wd, pmp_D, "Donut", r, phi, z, w0, z0, pump_E);  % 775[nm], z0 = 0
Ipump = pump.intensityProfileBLDumped(z);

% Pump beam propogating inside the sample:


% Shining pump beam on the sample, causing e-h generation and diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates spacial and temporal FCC distribution

% Probe laser:
prb_wd = 5e-11;   % pulse of 50[ps]
prb_D = sp.D1;
probe_E0 = sqrt(sp.pulse_energies(5));  % E0=1
probe = laser(sp.wl1, prb_wd, prb_D, "Gauss", r, phi, z, w0, 0, probe_E0);  % 1550nm

% Shining the probe on the sample:
%prb_prf_after = probeAbsorptionBS(probe, pDiff, r, z, t);






% Calculating the change in indexes at z=0 according to Bennet-Soref:
% Intensity slice at z=0, t=0  (adjust dims if needed)
% Ip0 = Ipump(:,1,1);                 % Nr×1
% 
% % Convert intensity -> carrier density (Nr×1)
% Ne = intensityToCarriers(Ip0, sp.alpha775, pump.lambda);
% Nh = Ne;
% 
% % Bennett–Soref (should work elementwise on vectors)
% [dn, dalpha] = BennetSoref(Ne, Nh, pump.lambda);   % Nr×1 each
% 
% % r: Nr×1 from 0..rmax
% %r = r(:);
% 
% % Map: r = |x|  (slice at y=0)
% dn_x     = interp1(r, real(dn(:)),     abs(x), 'linear', 0);
% dalpha_x = interp1(r, dalpha(:), abs(x), 'linear', 0);
% 
% n_x     = sp.n        + dn_x;
% alpha_x = sp.alpha775 + dalpha_x;
% 
% figure;
% yyaxis left;  plot(x*1e6, n_x);     ylabel('n');
% yyaxis right; plot(x*1e6, alpha_x*1e-2); ylabel('\alpha [1/cm]');
% xlabel('x [\mum]');
% title('n(x) and \alpha(x) at y=0, z=0, t=0 for Gaussian Pump');

% Displaying the FCC distribution in x-z for different time (diffusion):
% FCCPlot_xz(pDiff, Nr, Nz, r, x, z);

%iz = 1;         % z(1) = 0
% it = 1;         % t(1) = 0
% fprintf('z(iz) = %d', z(iz));

% FCC distribution at different times:
% figure;
% subplot(1,3,1);
% [P, xPlot, yPlot] = cylToCart(pDiff(:,:,1), r, iz);
% imagesc(xPlot, yPlot, P);
% colorbar;
% subplot(1,3,2);
% [P, xPlot, yPlot] = cylToCart(pDiff(:,:,(Nt+1)/2), r, iz);
% imagesc(xPlot, yPlot, P);
% colorbar;
% subplot(1,3,3);
% [P, xPlot, yPlot] = cylToCart(pDiff(:,:,Nt), r, iz);
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

% [P1, xPlot1, yPlot1] = cylToCart(prb_prf_after(:,:,1), r, 1);
% I1 = (abs(P1).^2) .* exp(-z(1)/sp.penDepth);
% [P2, xPlot2, yPlot2] = cylToCart(prb_prf_after(:,:,(Nt+1)/2), r, 1);
% I2 = (abs(P2).^2) .* exp(-z(1)/sp.penDepth);
% figure;
% subplot(1,2,1);
% surf(xPlot1*1e6, yPlot1*1e6, I1);
% shading interp
% xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [au]')
% colorbar
% subplot(1,2,2);
% surf(xPlot2*1e6, yPlot2*1e6, I2);
% shading interp
% xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('Intensity [au]')
% colorbar

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

