%% FINISHED

function CheckZ0(pump, z, r, x, t, probe, z0_vec,l)
    % Checks z0's effect on FWHM
    % ---------------------------------------------------------------------
    % This function changes z0 of the pump and checks the effect it has on
    % the probe's FWHM at the end of the sample, in two cases for the
    % probe's z0
    %
    % Plots the results
    % =====================================================================
    % INPUTS:
    %        pump - pump laser beam, Laser-type object
    %        z - z coordinate, propagation vector [m]
    %        r - radial spatial coordinate vector [m] 
    %        x - x coordinate vector [m]
    %        t - time vector [s]
    %        probe - probe laser beam, Laser-type object
    %        z0_vec - vector of z0 to check
    %        l - LG polynomial index
    % *********************************************************************
    
    sp = systemParameters();
    
    Nz = numel(z);
    %Nr = numel(r);
    Nx = numel(x);
    Nt = numel(t);
    Nz0 = numel(z0_vec);
    
    phi = atan(1);

    Ipump1_xz = zeros(Nx,Nz);
    %Iprobe1_xz = zeros(Nx,Nz,Nz0); Iprobe2_xz = zeros(Nx,Nz,Nz0);
    Ipropagate1_xz = zeros(Nx,Nz,Nt); %Ipropagate2_xz = zeros(Nx,Nz,Nt);
    %n_complex1 = zeros(Nr,Nz,Nt,Nz0); n_complex2 = zeros(Nr,Nz,Nt,Nz0);

    fwhm_end = zeros(2,Nz0); fwhm_end_pump = zeros(Nz0);
    
    % PML mask parameters:
    a = 4;
    b = 8;
    
    for i = 1:Nz0
        % Counter:
        disp(i);
        
        % Create new pump with the new z0:
        pump_i = laser(pump.lambda, pump.pulse_width, pump.pulse_energy, "Donut", r, phi, z, pump.w0, z0_vec(i),sp.n,l);
        Ipump1_xz = cylToCart(pump_i.intensityProfileBLDumped(z),r,x);
        %PF_plot_xz(Ipump1_xz*1e-4, z, x, sprintf('pump z0 = %.1f \\mum',z0_vec(i)*1e6));
        [fwhm_end_pump(i),~] = PF_x(Ipump1_xz(:,end)*1e-4,x,sprintf('Back Surafce: pump z0 = %.1f \\mum',z0_vec(i)*1e6),'Intensity [W/cm^2]');

        % Probe with constant z0=0 (sample surface):
        probe1 = laser(probe.lambda, probe.pulse_width, probe.pulse_energy, "Gauss", r, phi, z, probe.w0, probe.z0,sp.n);
        %Iprobe1_xz = cylToCart(probe1.intensityProfileBLDumped(z),r,x);
        
        % Same z0 as the pump's:
        %probe2 = laser(probe.lambda, probe.pulse_width, probe.pulse_energy, "Gauss", r, phi, z, probe.w0,  z0_vec(i),sp.n);
        %Iprobe2_xz = cylToCart(probe2.intensityProfileBLDumped(z),r,x);

        [pDiff1,~,~,Ipropagate1_xz] = runSimulation(x,z,r,t,pump_i,probe1,a,b);
        %[pDiff2,~,~,Ipropagate2_xz] = runSimulation(x,z,r,t,pump_i,probe2,a,b);
        [~, ~, it1] = findMax(pDiff1);
        %[~, ~, it2] = findMax(pDiff2);
        
        title1 = sprintf('pump z0 = %.1f \\mum ; probe z0 = %.1f \\mum',z0_vec(i)*1e6, 0);
        %PF_x_alongz(Ipropagate1_xz(:,:,it)*1e-4, x, z,title1,'Intensity [W/cm^2]');
        [fwhm_end(1,i),~] = PF_x(Ipropagate1_xz(:,end,it1)*1e-4,x,"Propagated Probe Intensity at Sample End "+title1,'Intensity [W/cm^2]');
        
        %title2 = sprintf('pump z0 = %.1f \\mum ; probe z0 = %.1f \\mum',z0_vec(i)*1e6, z0_vec(i)*1e6);
        %PF_x_alongz(Ipropagate2_xz(:,:,it)*1e-4, x, z,title2,'Intensity [W/cm^2]');       
        %[fwhm_end(2,i),~] = PF_x(Ipropagate2_xz(:,end,it2)*1e-4,x,"Propagated Probe Intensity at Sample End "+title2,'Intensity [W/cm^2]');
    end

    % Checking the FWHM at the end:
    figure('Color','w');
    plot(z0_vec*1e6,fwhm_end_pump*1e6, LineWidth=2); hold on;
    yline(0.91, 'r--');                   % Experimental result
    xline(0, 'k--'); xline(20, 'k--'); hold off;
    axis tight; grid on;
    xlabel('z0 [\mum]'); ylabel('Pump FWHM [\mum]');
    title("Pump's FWHM at the End of the Sample vs. Pump Waist Location on z Axis (z0)");
    
    figure('Color','w');
    plot(z0_vec*1e6,fwhm_end(1,:)*1e6, 'b', LineWidth=2); hold on;
    %plot(z0_vec*1e6,fwhm_end(2,:)*1e6, 'r:', LineWidth=2); hold on;
    yline(0.96, 'r--');                   % Experimental result
    xline(0, 'k--'); xline(20, 'k--'); hold off;
    axis tight; grid on;
    xlabel('z0 [\mum]'); ylabel('Probe FWHM [\mum]');
    legend('z_0probe=0','z_0probe=z_0pump');
    title("Probe's FWHM at the End of the Sample vs. Pump Waist Location on z Axis (z0)");
    
    figure('Color','w');
    scatter(z0_vec*1e6,abs(fwhm_end(1,:)-fwhm_end(2,:))*1e6, 'm'); hold on;
    xline(0, 'k--'); xline(20, 'k--'); hold off;
    axis tight; grid on;
    xlabel('z0 [\mum]'); ylabel('FWHM Difference [\mum]');
    title("Difference Between the Probe's FWHM in the Two Cases");
end