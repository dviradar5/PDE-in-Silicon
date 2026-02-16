function PF_FCCPlot_xz(pDiff, Nt, Nz, r, x, z)
    % Displaying the FCC distribution in x-z for different time (diffusion)
    % *********************************************************************

    rVec = r(:);
    xVec = [-flipud(rVec(2:end)); rVec];
    absx = abs(xVec);
    
    % index map: for each |x| find its r index (exact match)
    [~, idx] = ismembertol(absx, rVec, 0);  % tolerance 0 assumes exact; use 1e-12 if needed
    
    pCart = zeros(numel(xVec), Nz, Nt);
    for tt = 1:Nt
        pr = pDiff(:,:,tt);     % Nr×Nz
        pCart(:,:,tt) = pr(idx, :);
    end
    
    p1 = pCart(:,:,1); p2 = pCart(:,:,20); p3 = pCart(:,:,Nt);
    
    figure;
    
    subplot(1,3,1);
    imagesc(z*1e6, x*1e6, p1);
    axis xy
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title("t = 0");
    colorbar;
    
    subplot(1,3,2);
    imagesc(z*1e6, x*1e6, p2);
    axis xy
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title("t = 100[ps]");
    colorbar;
    
    subplot(1,3,3);
    imagesc(z*1e6, x*1e6, p3);
    axis xy
    xlabel('z [\mum]'); ylabel('x [\mum]');
    title("t = 2[ns]");
    colorbar;
    
    sgtitle("FCC x-z Distribution After Donut Pump at Different Times");
end