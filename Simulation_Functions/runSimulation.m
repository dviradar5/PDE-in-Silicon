%% FINISHED

function [pDiff,n_complex,Ipropagate,Ipropagate_xz] = runSimulation(x,z,r,t,pump,probe)
    % Quick Simulation Run
    % ---------------------------------------------------------------------
    % This function runs a quick simulation, without any plots and other
    % calculations
    %
    % Returns the FCC, refractive index and propogated Probe profile
    % matrices
    % =====================================================================
    % INPUTS:
    %        x - x coordinate vector [m]
    %        z - z coordinate, propagation vector [m]
    %        r - radial spatial coordinate vector [m] 
    %        t - time vector [s]
    %        pump - pump laser beam, Laser-type object
    %        probe - probe laser beam, Laser-type object
    % OUTPUTS:
    %        pDiff - FCC distribution, Nr x Nz x Nt, [1/m^3]
    %        n_complex - complex refractive index matrix at specific time, Nr x Nz
    %        Ipropagate - intensity in the sample, Nr x Nz [W/m^2]
    %        Ipropagate_xz - intensity in the sample, Nx x Nz [W/m^2]
    % *********************************************************************

    Nx = numel(x);
    Nr = numel(r);
    Nz = numel(z);
    Nt = numel(t);
    
    % Shining pump beam on the sample, causing e-h generation and diffusion:
    pDiff = FCCDiffusion(pump, t, r, z);    % Creates FCC distribution p(r,z,t)
        
    % Calculating the complex refractive index n(r,z,t) (changes due to FCC generation):
    n_complex = complexRefractiveIndex(pDiff, pump.lambda);
        
    % Probe Propogation:
    Ipropagate = zeros(Nr,Nz,Nt);       % Probe intensity after propogation
    Ipropagate_xz = zeros(Nx,Nz,Nt);
    
    for i = 1:Nt
        [~, Ipropagate(:,:,i)] = propagationBPM_rz(probe.profile(:,1), r, z, probe, n_complex(:,:,i));
        Ipropagate_xz(:,:,i) = cylToCart(Ipropagate(:,:,i),r,x);
    end

end