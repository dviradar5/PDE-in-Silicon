%% FINISHED

function [ad_z,iz] = absorptionDepth(I, z, lambda)
    % Calculates the Absorption Depth
    % ---------------------------------------------------------------------
    % Calculates absorption depth and calculates the theoretical value and
    % error
    % =====================================================================
    % INPUTS:
    %        I - intensity vector at r=0, Nz, [W/m^2]
    %        z - z coordinate vector [m]
    % OUTPUTS:
    %        iz - absorption depth index
    %        ad_z - absorption depth, [m]
    % *********************************************************************
    
    sp = systemParameters();

    % Finds the first crossing the I0/e value:
    iz = find(I <= I(1)/exp(1),1,"first");
    ad_z = z(iz);
    
    ad_theory = 1/sp.alpha;

    fprintf("Absorption Depth: %.2f [micron]\n", ad_z*1e6)
    fprintf("Theoretical Absorption Depth at %.0f[nm]: %.2f[micron]\n", lambda*1e9, ad_theory*1e6)
    fprintf("Error: %.2f [micron]\n", abs(ad_z-ad_theory)*1e6)

end