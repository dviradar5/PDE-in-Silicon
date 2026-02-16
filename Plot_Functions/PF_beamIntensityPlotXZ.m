%% DOCUMENT

function PF_beamIntensityPlotXZ(Ipump, r, z, x)
    zVec = z(:);
    xVec = x(:);
    
    absx = abs(xVec);
    
    IxZ = zeros(numel(xVec), numel(zVec));
    for k = 1:numel(zVec)
        IxZ(:,k) = interp1(r(:), Ipump(:,k), absx, 'linear', 'extrap');
    end
    
    figure;
    imagesc(zVec*1e6, xVec*1e6, IxZ);
    axis xy;
    xlabel('z [mum]');
    ylabel('x [mum]');
    title('Normalized I(x,z)');
    colorbar;
end