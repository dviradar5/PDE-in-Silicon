%% FINISHED

function CheckLGOrder(pump,z,r,x,t,probe,l_vec)
    % Checks LG order effect on FWHM
    % ---------------------------------------------------------------------
    % This function changes the order of the Laguerre polynomial of the
    % pump and checks the effect it has on its FWHM at the end of the
    % sample
    % =====================================================================
    % INPUTS:
    %        pump - pump laser beam, Laser-type object
    %        z - z coordinate, propagation vector [m]
    %        r - radial spatial coordinate vector [m]
    %        x - x coordinate vector [m]
    %        t - time vector [s]
    %        probe - probe laser beam, Laser-type object
    %        l_vec - vector of LG orders to check
    % *********************************************************************
    
    phi = atan(1);
    % PML mask parameters:
    a = 4;
    b = 8;

    results = [];
    best.fwhm = inf;
    
    fwhm_pump_end = zeros(1,numel(l_vec));
    fwhm_probe_end = zeros(1,numel(l_vec));
    
    for il = 1:numel(l_vec)
    
        l_i = l_vec(il);
    
        pump_i = laser(pump.lambda, pump.pulse_width, pump.pulse_energy,"Donut", r, phi, z, pump.w0, pump.z0,l_i);
        Ipump_i_xz = cylToCart(pump_i.intensityProfileBLDumped(z),r,x);
    
        % Pump FWHM at end of sample:
        [fwhm_pump_end(il),~] = PF_x(Ipump_i_xz(:,end)*1e-4, x,"Pump at sample end", "Intensity [W/cm^2]");
    
        % Run simulation:
        [pDiff, ~, Ipropagate, ~] = runSimulation(x, z, r, t, pump_i, probe,a,b);
    
        % Best time index:
        [~, ~, itMax] = findMax(pDiff);
        Ipropagate_xz = cylToCart(Ipropagate(:,:,itMax),r,x);
    
        % Probe FWHM at end of sample:
        [fwhm_probe_end(il),~] = PF_x(Ipropagate_xz(:,end)*1e-4, x);
    
        results = [results; l_i, fwhm_pump_end(il), fwhm_probe_end(il)];
    
        fprintf("\nLG order (l) = %d | pump = %.3f um, probe = %.3f um", ...
            l_i, fwhm_pump_end(il)*1e6 , fwhm_probe_end(il)*1e6);
    
        if fwhm_probe_end(il) < best.fwhm
            best.l = l_i;
            best.fwhm = fwhm_probe_end(il);
            best.fwhm_pump = fwhm_pump_end(il);
            best.fwhm_probe = fwhm_probe_end(il);
        end
    end
    
    best.results = array2table(results,'VariableNames', {'order','fwhm_pump','fwhm_probe_m'});
    
    fprintf("\nBEST RESULT:")
    fprintf("\nLG order (l) = %d", best.l)
    fprintf("Pump FWHM = %.3f [um]\n", best.fwhm_pump*1e6)
    fprintf("Probe FWHM = %.3f [um]\n", best.fwhm_probe*1e6)
       
    figure('Color','w');
    plot(l_vec, fwhm_pump_end*1e6, 'm-o', 'LineWidth', 1.5); hold on;
    plot(l_vec, fwhm_probe_end*1e6, 'g-o', 'LineWidth', 1.5);
    plot(best.l, best.fwhm_pump*1e6, 'mp','MarkerSize', 12, 'MarkerFaceColor', 'm');
    plot(best.l, best.fwhm_probe*1e6, 'gp','MarkerSize', 12, 'MarkerFaceColor', 'g');
    hold off; grid on; axis tight; 
    legend("Pump","Probe","Best Pump","Best Probe");
    xlabel('LG order |l|'); ylabel('FWHM [\mum]');
    title("Pump & Propagated Probe FWHM at Sample End vs. Pump Vortex Order");
end