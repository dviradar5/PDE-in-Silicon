%% demo_compare_solvers.m
% Minimal test comparing BPM and angular-spectrum split-step.
% The FDFD solver is included but should be run only on small/coarse grids.

clear; close all; clc;

probe.lambda = 775e-9;    % [m]
probe.n0 = 3.68;          % approximate real Si index reference; set yours

Nr = 201;
Nz = 251;
r = linspace(0, 8e-6, Nr).';
z = linspace(0, 25e-6, Nz).';

w0 = 2.378e-6;
E0_r = exp(-(r/w0).^2);   % field amplitude, so intensity radius differs

% Example uniform complex index. Replace this with your PDE-updated n_complex.
n_si = 3.68 + 1i*0.006;
n_complex = n_si * ones(Nr, Nz);

% Optional radial boundary sponge. Increase a only if edge reflections appear.
a = 6;
b = 12;

[E_bpm, I_bpm] = propagateField_rz('bpm',  E0_r, r, z, probe, n_complex, a, b);
[E_ssf, I_ssf] = propagateField_rz('ssfm', E0_r, r, z, probe, n_complex, a, b);

figure('Color','w');
plot(r*1e6, I_bpm(:,end)/max(I_bpm(:,end)), 'LineWidth', 2); hold on;
plot(r*1e6, I_ssf(:,end)/max(I_ssf(:,end)), '--', 'LineWidth', 2);
xlabel('r [um]'); ylabel('normalized intensity'); grid on;
legend('BPM', 'angular-spectrum split-step');
title('Output radial profile');

figure('Color','w');
imagesc(z*1e6, r*1e6, I_bpm/max(I_bpm(:))); axis xy; colorbar;
xlabel('z [um]'); ylabel('r [um]'); title('BPM normalized intensity');

%% Optional FDFD sanity check on much smaller grid
% Nr2 = 61; Nz2 = 81;
% r2 = linspace(0, 5e-6, Nr2).';
% z2 = linspace(0, 10e-6, Nz2).';
% E02 = exp(-(r2/w0).^2);
% n2 = n_si*ones(Nr2,Nz2);
% [E_fdfd, I_fdfd, info_fdfd] = propagateField_rz('fdfd', E02, r2, z2, probe, n2, 2, 8);
% figure('Color','w'); imagesc(z2*1e6, r2*1e6, I_fdfd/max(I_fdfd(:))); axis xy; colorbar;
% xlabel('z [um]'); ylabel('r [um]'); title('Scalar Helmholtz FDFD intensity');
