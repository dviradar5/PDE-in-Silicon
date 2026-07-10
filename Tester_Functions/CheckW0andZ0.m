%% CHECKKKKKKKKKKKKKK

function CheckW0andZ0(pump, z, r, x, t, probe, z0_vec, w0_vec, fwhm_pump, fwhm_probe, fileName)
    % ---------------------------------------------------------------------
    % Evaluates the effect of pump waist (w0) and focal position (z0) on 
    % the final propagated probe FWHM, and writes the results to a .txt
    % file
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
    %        fileName - text file name, string
    % *********************************************************************
    
    if nargin < 11 || isempty(fileName)
        fileName = 'Figures_and_Results\sweep.txt';
    end

    % Opening a text file with a header:
    fileID = fopen(fileName,'w+');
    fprintf(fileID,"w0 [um], z0 [um] | pump FWHM [um] | probe FWHM [um] | error [um]\n");

    sp = systemParameters();

    phi = atan(1);
    
    % PML mask parameters:
    a = 4;
    b = 8;
    
    % Initialize matrices where rows = w0, columns = z0
    errorMap = nan(numel(w0_vec), numel(z0_vec));
    probeFwhmMap = nan(numel(w0_vec), numel(z0_vec));
    results = [];
    
    % Track overall best performance
    best.error = inf;

    for iw0 = 1:numel(w0_vec)
        for iz0 = 1:numel(z0_vec)
            w0_i = w0_vec(iw0);
            z0_i = z0_vec(iz0);
            
            % Generate Pump 
            pump_i = laser(pump.lambda, pump.pulse_width, pump.pulse_energy, "Donut", r, phi, z, w0_i, z0_i,sp.n,1);
            Ipump_i_xz = cylToCart(pump_i.intensityProfileBLDumped(z), r, x);
            
            % Pump FWHM at end of sample
            % (Assuming PF_x returns FWHM without halting the loop)
            [fwhm_pump_end, ~] = PF_x(Ipump_i_xz(:,end)*1e-4, x, "Pump at sample end", "Intensity [W/cm^2]");
            
            probe1 = laser(probe.lambda, probe.pulse_width, probe.pulse_energy,"Gauss", r, phi, z, probe.w0, probe.z0,sp.n);

            % Run simulation
            [pDiff, ~, Ipropagate, ~] = runSimulation(x, z, r, t, pump_i, probe1, a, b);
            
            % Find peak plasma time index
            [~, ~, itMax] = findMax(pDiff);
            
            % Probe FWHM at end of sample
            % Checking dimensions to prevent errors if Ipropagate is 2D vs 3D
            if ndims(Ipropagate) == 3
                probe_fwhm = FWHM(Ipropagate(:,:,itMax), r);
            else
                probe_fwhm = FWHM(Ipropagate, r);
            end
            fwhm_probe_end = probe_fwhm(end);
           
            % 2. Calculate objective error to minimize
            err_pump = abs(fwhm_pump_end - fwhm_pump);
            err_probe = abs(fwhm_probe_end - fwhm_probe);
            total_err = err_pump + err_probe;
            
            % Store values in the 2D matrices
            errorMap(iw0, iz0) = total_err;
            probeFwhmMap(iw0, iz0) = fwhm_probe_end;
            
            results = [results; w0_i, z0_i, fwhm_pump_end, fwhm_probe_end, total_err]; 
            
            fprintf(fileID,"%.3f, %.1f | %.3f | %.3f | %.3f\n", ...
                w0_i*1e6, z0_i*1e6, fwhm_pump_end*1e6, fwhm_probe_end*1e6, total_err*1e6);
            
            % Track global minimum
            if total_err < best.error
                best.w0 = w0_i;
                best.z0 = z0_i;
                best.iz0 = iz0; % Store the matrix column index of the best z0
                best.fwhm_pump = fwhm_pump_end;
                best.fwhm_probe = fwhm_probe_end;
                best.error = total_err;
            end
        end
    end
    
    % Store and display results table
    best.results = array2table(results, 'VariableNames', {'w0', 'z0', 'fwhm_pump', 'fwhm_probe', 'error'});
    fprintf(fileID,"\nBEST OVERALL RESULT:\n");
    fprintf(fileID,"w0 = %.3f [um]\n", best.w0*1e6);
    fprintf(fileID,"z0 = %.3f [um]\n", best.z0*1e6);
    fprintf(fileID,"Pump FWHM = %.3f [um]\n", best.fwhm_pump*1e6);
    fprintf(fileID,"Probe FWHM = %.3f [um]\n", best.fwhm_probe*1e6);
    fprintf(fileID,"Error = %.3f [um]\n", best.error*1e6);
    
    % Closing the file:
    fclose(fileID);

    % =====================================================================
    % Data Extraction for FWHM Plots
    % =====================================================================
    
    % Curve 1: Probe FWHM vs w0 for the SINGLE best global z0
    % We extract the entire column corresponding to 'best.iz0'
    probe_fwhm_fixed_z0 = probeFwhmMap(:, best.iz0);
    
    % Curve 2: Probe FWHM vs w0, optimizing z0 independently for each w0
    % min(matrix, [], 2) finds the minimum across columns (the z0 dimension)
    [~, best_iz0_per_w0] = min(errorMap, [], 2); 
    probe_fwhm_optimal_z0 = zeros(numel(w0_vec), 1);
    
    for iw0 = 1:numel(w0_vec)
        % For each w0 row, pull the FWHM at its specific best z0 column
        probe_fwhm_optimal_z0(iw0) = probeFwhmMap(iw0, best_iz0_per_w0(iw0));
    end

    % =====================================================================
    % Plotting
    % =====================================================================
    
    % Plot A: Error Colormap (Your original requested plot)
    figure('Color','w');
    imagesc(z0_vec*1e6, w0_vec*1e6, errorMap*1e6); hold on;
    plot(best.z0*1e6, best.w0*1e6, 'rx', 'DisplayName', "Minimum", 'MarkerSize', 15, 'LineWidth', 2); hold off;
    axis xy; colorbar;
    xlabel('z_0 [\mum]'); ylabel('w_0 [\mum]');
    title("Total Error vs. z_0 and w_0");    
    
    % Plot B: Probe FWHM Comparisons
    figure('Color', 'w');
    plot(w0_vec*1e6, probe_fwhm_fixed_z0*1e6, '-o', 'LineWidth', 2, ...
        'DisplayName', sprintf('Fixed optimal z_0 = %.1f \\mum', best.z0*1e6));
    hold on;
    plot(w0_vec*1e6, probe_fwhm_optimal_z0*1e6, '-s', 'LineWidth', 2, ...
        'DisplayName', 'Independently optimized z_0 per w_0');
    hold off;
    
    grid on;
    xlabel('Pump Waist w_0 [\mum]');
    ylabel('Propagated Probe FWHM [\mum]');
    title('Probe Focus Behavior vs. Pump Waist');
    legend('Location', 'best');
end