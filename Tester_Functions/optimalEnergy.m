%% FINISHED

function [opt_energy,opt_fwhm] = optimalEnergy(pump,z,r,x,t,probe)
    % Checks for the Best Fitting Energy
    % ---------------------------------------------------------------------
    % This function changes the pump's energy in search for the erengy that
    % would bring the focal point (minimal FWHM) of the probe as close as 
    % possible to the back end of the sample
    % =====================================================================
    % INPUTS:
    %        pump - pump laser beam, Laser-type object
    %        z - z coordinate, propagation vector [m]
    %        r - radial spatial coordinate vector [m] 
    %        x - x coordinate vector [m]
    %        t - time vector [s]
    %        probe - probe laser beam, Laser-type object
    % OUTPUTS:
    %        opt_energy - the best fitting energy, [nJ]
    %        opt_fwhm - FWHM at the sample end, [micron]
    % *********************************************************************

    energies = (1:1:100)*1e-9;       % Change to fit destruction limit
    opt_energy = energies(1);
    opt_fwhm = 0;
    opt_z = 0;
    
    % PML mask parameters:
    a = 4;
    b = 8;

    maxI_vec = zeros(1,numel(energies));
    maxILoc_vec = zeros(1,numel(energies));

    FWHMLoc_vec = zeros(1,numel(energies));
    minFWHM_vec = zeros(1,numel(energies));
    endFWHM_vec = zeros(1,numel(energies));
    
    bestDist = inf;
    
    for i = 1:numel(energies) 
        % Create new pump and probe:
        pump1 = laser(pump.lambda, pump.pulse_width, energies(i), "Donut", r, atan(1), z, pump.w0, pump.z0,l);
        probe1 = laser(probe.lambda, probe.pulse_width, energies(i)/100, "Gauss", r, atan(1), z, probe.w0, 0);
        
        [~,~,I,~] = runSimulation(x,z,r,t,pump1,probe1,a,b);
        Irz = I(:,:,11);

        fwhm_z = FWHM(Irz, r);
        
        [minFWHM, j] = min(fwhm_z);
        FWHMLoc_vec(i) = z(j);
        minFWHM_vec(i) = minFWHM;
        endFWHM_vec(i) = fwhm_z(end);
        
        [maxI_vec(i), linIdx] = max(Irz(:));
        [~, kZ] = ind2sub(size(Irz), linIdx);
        maxILoc_vec(i) = z(kZ);

        if abs(z(j)-z(end)) < bestDist
            bestDist = abs(z(j)-z(end));
            opt_energy = energies(i);
            opt_fwhm = minFWHM;
            opt_z = z(j);
        end
    end

    fprintf("\nOptimal Energy: %.2f[nJ], at z = %.2f[um]", opt_energy*1e9, opt_z*1e6);
    fprintf("\nFWHM: %.2f[um]",opt_fwhm*1e6);
    
    figure('Color','w');
    plot(energies*1e9, FWHMLoc_vec*1e6, '-o', 'LineWidth', 1.5);
    grid on; axis tight;
    xlabel('Pump Energy [nJ]'); ylabel('zFWHM [\mum]');
    title('Minimal Probe FWHM z-axis locztion vs. Pump Energy');

    figure('Color','w');
    plot(energies*1e9, minFWHM_vec*1e6, '-o', 'LineWidth', 1.5);
    grid on; axis tight;
    xlabel('Pump Energy [nJ]'); ylabel('Minimal FWHM [\mum]');
    title('Minimal Probe FWHM vs. Pump Energy');

    figure('Color','w');
    plot(energies*1e9, endFWHM_vec*1e6, '-o', 'LineWidth', 1.5);
    grid on; axis tight;
    xlabel('Pump Energy [nJ]'); ylabel('FWHM at Sample End [\mum]');
    title('Probe FWHM at Sample End vs. Pump Energy');
        
    figure('Color','w');
    plot(energies*1e9, maxI_vec*1e-4, '-o', 'LineWidth', 1.5);
    grid on; axis tight;
    xlabel('Pump Energy [nJ]'); ylabel('Max Intensity [W/cm^2]');
    title('Probe Max Intensity vs. Pump Energy');

    figure('Color','w');
    plot(energies*1e9, maxILoc_vec*1e6, '-o', 'LineWidth', 1.5);
    grid on; axis tight;
    xlabel('Pump Energy [nJ]'); ylabel('zMax [\mum]');
    title('Probe Max Intensity Location vs. Pump Energy');
end