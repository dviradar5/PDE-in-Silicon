%% FIXXXXXXXXXXXXXXX

function [deltaFWHM_all, minFWHM_all, izSlices] = sliceContribution(pDiff, pump, probe, r, z, width)
    % Slicing Function
    % ---------------------------------------------------------------------
    % This function creates slices of the sample in order to check each
    % one's contribution to the overall FWHM
    % =====================================================================
    % INPUTS:
    %        pDiff - FCC distribution, Nr x Nz x Nt, [1/m^3]
    %        pump - pump laser beam, Laser-type object
    %        probe - probe laser beam, Laser-type object
    %        z - z coordinate, propagation vector [m]
    %        r - radial spatial coordinate vector [m] 
    %        width - each slices' width [microns]
    % OUTPUTS:
    %        deltaFWHM_all - FWHM distance vector from the original FWHM, size: number of slices
    %        minFWHM_all - all the slices' minimal FWHM vector, size: number of slices
    %        izSlices - slices' indecies vector, size: number of slices
    % *********************************************************************
    sp = systemParameters();

    Nr = size(pDiff,1);
    Nz = size(pDiff,2);
    
    % Translating width to indices:
    if (sp.Lz*1e6/Nz > width)
        error("Width too small; needs to be bigger than the grid resolution, %f micron", 1e6*sp.Lz/Nz);
    end
    
    dzStep = max(1, round(width*Nz/(sp.Lz*1e6)));
    
    if dzStep >= Nz
        izSlices = 1;       % 1 slice case (regular calculation)
    else
        izSlices = 1:dzStep:Nz;
    end
    
    n_base = complexRefractiveIndex(sp.ni * ones(Nr,Nz), pump.lambda);
    [~, I_base] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_base);
    
    FWHM_base = FWHM(I_base, r);
    minFWHM_base = min(FWHM_base);
    
    Ns = numel(izSlices);
    minFWHM_all = zeros(1,Ns);
    zmin = zeros(1,Ns);
    deltaFWHM_all = zeros(1,Ns);
    
    for k = 1:Ns

        iz = izSlices(k);
    
        if dzStep >= Nz
            s = pDiff;
        else
            s = sp.ni * ones(Nr,Nz);        
            
            % Each slice has a width:
            for l = 0:dzStep-1
                if (iz+l <= Nz)
                    s(:,iz+l) = pDiff(:,iz+l);
                end
            end
        end
        
    
        n_slice = complexRefractiveIndex(s, pump.lambda);
        [~, I_slice] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_slice);
    
        FWHM_slice = FWHM(I_slice, r);          % Nz vector

        [minFWHM_all(k),i] = min(FWHM_slice(:));
        zmin(k) = z(i);

        deltaFWHM_all(k) = minFWHM_base - minFWHM_all(k);
    end
    
    % Results:
    [bestDelta, bestIdx] = max(deltaFWHM_all);
    best_iz = izSlices(bestIdx);
    
    fprintf("Number of Slices: %d\n", Ns)
    fprintf("Slice width: %.2f micron\n", width)
    fprintf('Best contributing slice: #%d, z = %.2f um\n', bestIdx, z(best_iz)*1e6);
    fprintf('Minimum FWHM: %.4e m, Achieved at z = %.2f um\n', minFWHM_all(bestIdx), zmin(bestIdx)*1e6);
    fprintf('Improvement: %.4e m\n', bestDelta);

    % Plotting FWHM:
    figure;
    
    subplot(2,1,1);
    scatter(1:Ns,minFWHM_all*1e6);
    hold on;
    yline(minFWHM_base*1e6,'r--');
    hold off;
    axis tight; grid on;
    ylabel('FWHM [\mum]');
    legend("Slice's Min. FWHM","Original Min. FWHM");
    title('Minimal FWHM and Location (Focal Length)');
    
    subplot(2,1,2);
    scatter(1:Ns,zmin*1e6);
    hold on;
    yline(z(end)*1e6,'r--');
    hold off;
    axis tight; grid on;
    legend("FWHM Location","Sample End");
    xlabel('Slices'); ylabel('Minimal FWHM Location [\mum]');
end