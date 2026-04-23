%% FINIDSHED

function [xMax,zMax,tMax] = findMax(A)
    % Maximum  value finder function
    % ---------------------------------------------------------------------
    % Finds the maximal value of the matrix and returns the indices
    % =====================================================================
    % INPUTS:
    %        A - matrix, Nx x Nz x Nt
    % OUTPUTS:
    %        xMax - 1st dimension index
    %        zMax - 2nd dimension index
    %        tMax - 3rd dimension index
    % *********************************************************************
    
    [~, B] = max(A(:));
    [~, ~, tMax] = ind2sub(size(A), B);
    
    slice = A(:,:,tMax);
    [~, B] = max(slice(:));
    [xMax, zMax] = ind2sub(size(slice), B);
end