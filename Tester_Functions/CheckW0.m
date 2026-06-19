%% FINISHED

function CheckW0(pump, z, r, x, t, probe, w0_vec, l)
    % Checks pump w0 effect on pump/probe FWHM at sample end
    % ---------------------------------------------------------------------
    % This function changes the waist radius w0 of the pump and checks the
    % effect it has on the probe's FWHM at the end of the sample,
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
    %        w0_vec - vector of w0 to check
    %        l - LG polynomial index
    % *********************************************************************

    Nw0 = numel(w0_vec);

    phi = atan(1);

    fwhm_end_pump = zeros(1, Nw0);
    fwhm_end_probe = zeros(1, Nw0);
    
    inner_D = zeros(1, Nw0);
    outer_D = zeros(1, Nw0);

    % PML mask parameters:
    a = 4;
    b = 8;

    for iw = 1:Nw0
        % Create new pump with new w0:
        pump_i = laser(pump.lambda, pump.pulse_width, pump.pulse_energy,"Donut", r, phi, z, w0_vec(iw), pump.z0, l);
        Ipump = pump_i.intensityProfileBLDumped(z);
        Ipump_xz = cylToCart(Ipump, r, x);

        % Pump FWHM at sample end:
        [fwhm_end_pump(iw), ~] = PF_x(Ipump_xz(:,end)*1e-4, x, ...
            sprintf('Back Surface: pump w0 = %.3f \\mum', w0_vec(iw)*1e6),'Intensity [W/cm^2]');

        % Create probe:
        probe1 = laser(probe.lambda, probe.pulse_width, probe.pulse_energy,"Gauss", r, phi, z, probe.w0, 0);

        % Run simulation:
        [pDiff, ~, ~, Ipropagate_xz] = runSimulation(x, z, r, t, pump_i, probe1, a, b);

        % Find strongest focusing time:
        [~, ~, it] = findMax(pDiff);

        % Probe FWHM at sample end:
        [fwhm_end_probe(iw), ~] = PF_x(Ipropagate_xz(:,end,it)*1e-4, x, ...
            "Propagated Probe Intensity at Sample End ",'Intensity [W/cm^2]');
        
        % Pump diameters at the end:
        [~,r1,r2] = computeDiameter(Ipump(:,end), r, z(end), "Donut",l, pump.lambda, w0_vec(iw), pump.z0);
        inner_D(iw) = 2*r1(1);
        outer_D(iw) = 2*r2(1);
    end

    % Plots:
    % Pump diameters vs w0:
    figure('Color','w');
    subplot(1,2,1);
    plot(w0_vec*1e6, inner_D*1e6, 'bo-', 'LineWidth', 2);
    grid on; axis tight;
    xlabel('Pump w_0 [\mum]');
    ylabel('Inner Diameter [\mum]');
    
    subplot(1,2,2);
    plot(w0_vec*1e6, outer_D*1e6, 'ro-', 'LineWidth', 2);
    grid on; axis tight;
    xlabel('Pump w_0 [\mum]');
    ylabel('Outer Diameter [\mum]');
    
    sgtitle('Pump Diameters at Sample End vs. Pump Waist w_0');

    % Probe FWHM vs pump FWHM:
    figure('Color','w');
    plot(fwhm_end_pump*1e6, fwhm_end_probe*1e6, 'bo-', LineWidth=2);
    axis tight; grid on;
    xlabel('Pump FWHM [\mum]');
    ylabel('Probe FWHM [\mum]');
    title("Probe FWHM at Sample End vs. Pump FWHM for Different Pump {\omega}_0 Values");

    % FWHM vs w0:
    figure('Color','w');
    plot(w0_vec*1e6, fwhm_end_pump*1e6, 'k-', LineWidth=2); hold on;
    plot(w0_vec*1e6, fwhm_end_probe*1e6, 'b-', LineWidth=2);hold off;
    axis tight; grid on;
    xlabel('Pump w_0 [\mum]'); ylabel('FWHM [\mum]');
    legend('Pump', 'Probe', Location='best');
    title("FWHM at Sample End vs. Pump Waist w_0");
end