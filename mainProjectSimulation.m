%% DOCUMENT
clc; clear;
%close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp = systemParameters();

% Temporal vector intialization:
tf = 30000;                 % 30[ns]
t = 0:10:tf;                % Time vector [ps]
t = t .* 1e-12;             % Time vector, 0-30[ns], [s]
dt = t(2)-t(1);
Nt = numel(t);

% Spacial vector intialization:
Nr = 201;                               % Number of elements
Nz = 201;

z = linspace(0, sp.Lz, Nz);
r = linspace(0,sqrt(2)*sp.Lx,Nr);       % 10x10 micron sample so r covers it
phi = atan(1);                          % y = x, radial symmetry

x = linspace(-max(r), max(r), 2*Nr-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ENERGIES
% Laser constants:
pmp_wd = 3e-11;         % Pulse of 30[ps]
pmp_D = sp.D2;
w0 = 5e-6;
z0 = 0;                 % Beam waist location (-7 in old simulation)
pump_E = 16e-9;         % [J] FIXXXXXXXXXXXXX
%dE = 8e-9;             % Size of each energy portion [J]
%nsim = round(pump_E/dE);    % Number of iterations

% Pump laser:
pump = laser(sp.wl2, pmp_wd, pump_E, pmp_D, "Donut", r, phi, z, w0, z0);  % 775[nm]
Ipump = pump.intensityProfileBLDumped(z);
%beamIntensityPlotXZ(Ipump, r, z, x)
%beamPlotXY(r, x, Ipump(:,1), 1, 10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shining pump beam on the sample, causing e-h generation and diffusion:
pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution
%tIdx = [1, 20, 200, Nt];
%FCCPlot_xz_from_rzt(pDiff, r, z, t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAYBE MOVE TO FUNCTION
% Calculating the complex refractive index n(r,z,t) (it changes due to diffusion):
k0 = sp.alpha * pump.lambda / (4*pi);
n_complex = zeros(Nr, Nz, Nt);      % complex refractive index over time

for it = 1:Nt
    N = pDiff(:,:,it);              % Ne = Nh since 1 photon = e+h

    % Bennett–Soref expects Ne, Nh at pump wavelength:
    [dn, dalpha] = BennettSoref(N, N, pump.lambda);

    n_complex(:,:,it) = (sp.n + dn) + 1i*(k0 + dalpha * pump.lambda/ (4*pi));
    
    % Calculating relative permitivity for maxwell's equations later:
    %epsilon_r = n_complex^2;
end

%PF_complexRefractiveIndex(n_complex, r, z, x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probe laser:
prb_wd = 5e-11;   % pulse of 50[ps]
prb_D = sp.D1;
probe_E = sp.pulse_energies(5);         % [J]
probe = laser(sp.wl2, prb_wd, probe_E, prb_D, "Gauss", r, phi, z, w0, 0);   % 775[nm]

Ein_r = probe.profile(:,1);     % field at z=0 from your laser class
Iin_r = abs(Ein_r).^2;
Iin_r = Iin_r / max(Iin_r);

figure;
plot(r*1e6, Iin_r, 'LineWidth', 1.5);
xlabel('r [\mum]'); ylabel('norm |E|^2');
title('Input probe intensity at z=0 (should be Gaussian in r)');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- run maxwell for selected times only (faster to inspect) ---
t_idx = unique([1, round(Nt/4), round(Nt/2), Nt]);

for kk = 1:numel(t_idx)
    i = t_idx(kk);

    [E_rz, I_rz] = maxwell(r, z, probe.lambda, n_complex(:,:,i), probe.profile(:,1));
    %[E_rz, I_rz] = maxwell_fdfd(r, z, probe.lambda, n_complex(:,:,i), probe.profile(:,1));


    % Normalize globally for each snapshot
    Iplot = I_rz;% / max(I_rz(:));

    figure;
    imagesc(z*1e6, r*1e6, Iplot);
    set(gca,'YDir','normal'); axis tight;
    xlabel('z [\mum]'); ylabel('r [\mum]');
    title(sprintf('Probe intensity |E|^2 (normalized), t = %.3f ns', t(i)*1e9));
    colorbar;

    % Also show the surface (z=0) in Cartesian
    Ir0 = I_rz(:,1);
    %Ir0 = Ir0 / max(Ir0);

    % robust Cartesian mapping (avoid cylToCart if you want)
    x = linspace(-max(r), max(r), 2*numel(r)-1);
    [X,Y] = meshgrid(x,x);
    R = hypot(X,Y);
    Ixy = interp1(r, Ir0, R, 'linear', 0);

    figure;
    imagesc(x*1e6, x*1e6, Ixy);
    axis image; set(gca,'YDir','normal');
    xlabel('x [\mum]'); ylabel('y [\mum]');
    title(sprintf('Surface intensity (Cartesian), t = %.3f ns', t(i)*1e9));
    colorbar;

    % And a 1D cut (x-axis) to verify Gaussian shape
    mid = ceil(size(Ixy,1)/2);
    figure;
    plot(x*1e6, Ixy(mid,:), 'LineWidth', 1.5);
    xlabel('x [\mum]'); ylabel('norm I(x, y=0)');
    title(sprintf('Surface linecut (should be Gaussian), t = %.3f ns', t(i)*1e9));
    grid on;
end


probeComparison();



% Probe beam propogating inside the sample:
E_out = zeros(Nr,Nz,Nt);
I_out = zeros(Nr,Nz,Nt);
for i=1:Nt
    % solving maxwell in r-z for specific t:
    [E_out(:,:,i), I_out(:,:,i)] = maxwell(r, z, probe.lambda, n_complex(:,:,i), probe.profile(:,1));
end



%figure; plot(r*1e6, I_out(:,1,1)/max(I_out(:,1,1)));
%title('I(r) at z=0'); xlabel('r [\mum]'); ylabel('normalized I');

I1 = I_out(:,:,1);


figure;
imagesc(z, r, I1);
axis image; set(gca,'YDir','normal');
xlabel('z [\mum]'); ylabel('r [\mum]');
title('probe intensity, t=0');
colorbar;





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

