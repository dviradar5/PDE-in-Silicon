function FCCPlot_xz(pDiff, Nr, Nz, r, x, z)
    % Displaying the FCC distribution in x-z for different time (diffusion)
    % *********************************************************************

    p1 = zeros(2*Nr-1, Nz);             % pDiff(x,z), cartesian not cylindrical                          
    p2 = zeros(2*Nr-1, Nz);                     
    p3 = zeros(2*Nr-1, Nz);                         
    
    for iz = 1:Nz
        [P1,~,~] = cylToCart(pDiff(:,:,1), r, iz);      % t = 0
        [P2,~,~] = cylToCart(pDiff(:,:,11), r, iz);     % t = 100[ps]
        [P3,~,~] = cylToCart(pDiff(:,:,201), r, iz);    % t = 2[ns]
        p1(:,iz) = P1(:, Nr);      % y=0 at Nr in a (2*Nr-1) x (2*Nr-1) matrix
        p2(:,iz) = P2(:, Nr);
        p3(:,iz) = P3(:, Nr);
    end
    
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