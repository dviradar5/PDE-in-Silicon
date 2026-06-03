%% FINISHED

function CheckW0andZ0(pump, z, r, x, t, probe,z0_vec,w0_vec,fwhm_pump,fwhm_probe)
    % Checks w0's effect on FWHM
    % ---------------------------------------------------------------------
    % This function changes w0 of the pump and checks the effect it has on
    % the probe's FWHM at the end of the sample
    %
    % We also use CheckZ0 to see the effect they have together
    % =====================================================================
    % INPUTS:
    %        pump - pump laser beam, Laser-type object
    %        z - z coordinate, propagation vector [m]
    %        r - radial spatial coordinate vector [m] 
    %        x - x coordinate vector [m]
    %        t - time vector [s]
    %        probe - probe laser beam, Laser-type object
    %        z0_vec - vector of z0 to check
    %        w0_vec - vector of w0 to check
    %        fwhm_pump - desired pump fwhm at the sample end
    %        fwhm_probe - desired probe fwhm at the sample end
    % *********************************************************************
    
    phi = atan(1);

    results = [];

    best.error = inf;

    for iw0 = 1:numel(w0_vec)
        for iz0 = 1:numel(z0_vec)

            w0_i = w0_vec(iw0);
            z0_i = z0_vec(iz0);

            pump_i = laser(pump.lambda, pump.pulse_width, pump.pulse_energy,"Donut", r, phi, z, w0_i, z0_i);
            Ipump_i_xz = cylToCart(pump_i.intensityProfileBLDumped(z),r,x);

            % Run simulation:
            [pDiff, ~, Ipropagate, ~] = runSimulation(x, z, r, t, pump_i, probe);

            % Best time index:
            [~, ~, itMax] = findMax(pDiff);

            % Pump FWHM at end of sample:
            %pump_fwhm = FWHM(Ipump_i,r);
            %fwhm_pump_end = pump_fwhm(end);
            fwhm_pump_end = PF_x(Ipump_i_xz(:,end)*1e-4, x,"Pump at sample end", "Intensity [W/cm^2]");

            % Probe FWHM at end of sample:
            probe_fwhm = FWHM(Ipropagate(:,:,itMax),r);
            fwhm_probe_end = probe_fwhm(end);
           
            % Error in microns:
            err_pump = abs(fwhm_pump_end - fwhm_pump);
            err_probe = abs(fwhm_probe_end - fwhm_probe);

            total_err = err_pump + err_probe;
            errorMap(iz0,iw0) = total_err;
            
            results = [results; w0_i, z0_i, fwhm_pump_end, fwhm_probe_end, total_err];

            fprintf("w0 = %.3f um, z0 = %.1f um | pump = %.3f um, probe = %.3f um | error = %.3f um\n", ...
                w0_i*1e6, z0_i*1e6,fwhm_pump_end*1e6, fwhm_probe_end*1e6, total_err*1e6);

            if total_err < best.error
                best.w0 = w0_i;
                best.z0 = z0_i;
                best.fwhm_pump = fwhm_pump_end;
                best.fwhm_probe = fwhm_probe_end;
                best.error = total_err;
                best.itMax = itMax;
            end
        end
    end

    best.results = array2table(results,'VariableNames', {'w0_m','z0_m','fwhm_pump_m','fwhm_probe_m','error_m'});

    fprintf("\nBEST RESULT:\n")
    fprintf("w0 = %.3f [um]\n", best.w0*1e6)
    fprintf("z0 = %.3f [um]\n", best.z0*1e6)
    fprintf("Pump FWHM = %.3f [um]\n", best.fwhm_pump*1e6)
    fprintf("Probe FWHM = %.3f [um]\n", best.fwhm_probe*1e6)
    fprintf("Error = %.3f [um]\n", best.error*1e6)
    
    figure;
    imagesc(z0_vec*1e6,w0_vec*1e6,errorMap*1e6); hold on;
    plot(best.z0*1e6, best.w0*1e6,'x',DisplayName="Minimum", MarkerSize=40); hold off;
    axis xy; colorbar;
    xlabel('z0 [\mum]'); ylabel('w0 [\mum]');
    title("Total Error (Pump+Probe) vs. z0 and w0");

end