%% FINISHED

function [P_lost, T, P_back] = powerLoss(I, r)
    % Calculates Power Loss
    % ---------------------------------------------------------------------
    % Calculates the power difference between the first and last z values
    % of a given intensity matrix
    % =====================================================================
    % INPUTS:
    %        I - intensity matrix, Nr x Nz, [W/m^2]
    %        r - radial coordinate vector [m]
    % OUTPUT:
    %        P_lost - Pi-Pf, [W]
    %        T - percentage of power transmitted
    %        P_back - Pf, [W]
    % *********************************************************************

    I_front = I(:, 1);
    I_back = I(:, end);
    
    % Integrating over the whole r slice (including Jacobian):
    P_front = 2*pi*trapz(r(:), I_front.* r(:));
    P_back = 2*pi*trapz(r(:), I_back.*r(:));
    
    P_lost = P_front - P_back;                  % [W]
    T = (P_back / P_front)*100;                 % Transmission, percentage
    
    fprintf("\nPower lost: %.3f[W]",P_lost);
    fprintf("\nTransmission: %.3f%%\n",T);
end