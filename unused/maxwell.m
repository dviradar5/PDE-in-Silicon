%% DOCUMENT
%% BASED ON....

function [E_rz, I_rz] = maxwell(r, z, lambda, n_complex_rz, E0_r)
    % maxwell  Scalar paraxial Maxwell solver (BPM) in cylindrical symmetry.
    %
    % Solves for the slowly-varying envelope A(r,z) assuming exp(i*k0*n_ref*z).
    % Includes diffraction + phase shift from dn(r,z) + absorption from imag(n).
    %
    % Inputs:
    %   r            radial grid [m], Nr x 1
    %   z            propagation grid [m], 1 x Nz (or Nz x 1)
    %   lambda       wavelength [m]
    %   n_complex_rz complex index on grid (Nr x Nz): n - i*k
    %   E0_r         complex field at z=0 (Nr x 1)
    %
    % Outputs:
    %   E_rz         complex field E(r,z) (Nr x Nz)
    %   I_rz         intensity proxy |E|^2 (Nr x Nz)
    arguments
        r (:,1) double {mustBeNonnegative}
        z (:,1) double
        lambda (1,1) double {mustBePositive}
        n_complex_rz (:,:) double
        E0_r (:,1) double
    end
    
    sp=systemParameters();
    
    Nr = numel(r);
    Nz = numel(z);
    
    % Checking input validity:
    if size(n_complex_rz,1) ~= Nr || size(n_complex_rz,2) ~= Nz
        error("n_complex_rz must be size Nr x Nz to match r and z.");
    end
    if numel(E0_r) ~= Nr
        error("E0_r must be Nr x 1 to match r.");
    end

    dz = z(2) - z(1);
    
    if any(abs(diff(z) - dz) > 1e-15)
        error("z must be uniformly spaced for this BPM implementation.");
    end

    % Reference index for carrier/phase splitting
    n_ref = sp.n;
    k0    = 2*pi/lambda;

    % Build cylindrical transverse Laplacian operator:
    % (1/r) d/dr (r dA/dr)
    dr = r(2) - r(1);
    
    if any(abs(diff(r) - dr) > 1e-15)
        error("r must be uniformly spaced for this BPM implementation.");
    end

    main = zeros(Nr,1);
    up   = zeros(Nr,1);
    down = zeros(Nr,1);

    % Neumann at r=0: dA/dr = 0
    % Use symmetric FD at r=0
    main(1) = -2/dr^2;
    up(1)   =  2/dr^2;

    % Interior points
    for i = 2:Nr-1
        ri = r(i);
        main(i) = -2/dr^2;
        up(i)   =  (1/dr^2) + (1/(2*dr*ri));
        down(i) =  (1/dr^2) - (1/(2*dr*ri));
    end

    % Outer boundary at r=Rmax: simple absorbing-ish boundary
    % (Dirichlet A=0 is crude but stable)
    main(end) = 1;
    down(end) = 0;
    up(end)   = 0;

    L = spdiags([down main up], [-1 0 1], Nr, Nr);

    % Crank–Nicolson diffraction step matrices:
    c = 1i * dz / (4*k0*n_ref);
    A_mat = speye(Nr) - c*L;
    B_mat = speye(Nr) + c*L;

    % Optional soft absorber near boundary to reduce reflections:
    Rmax = r(end);
    absorber = ones(Nr,1);
    edge0 = 0.85*Rmax;
    idx = r > edge0;
    absorber(idx) = exp(-((r(idx)-edge0)/(Rmax-edge0)).^2 * 4);

    % Allocate outputs
    E_rz = zeros(Nr, Nz);
    E_rz(:,1) = E0_r;

    % Propagate
    for iz = 1:Nz-1
        n_here = n_complex_rz(:,iz);

        % Split-step "potential": phase & absorption relative to n_ref
        dn_eff = n_here - n_ref; % complex
        pot_half = exp(1i*k0*dn_eff*dz/2); % includes absorption if imag(dn_eff) != 0

        A0 = E_rz(:,iz);

        % Half potential
        A1 = A0 .* pot_half;

        % Diffraction CN step
        rhs = B_mat * A1;
        A2  = A_mat \ rhs;

        % Half potential again
        n_next = n_complex_rz(:,iz+1);
        dn_eff_next = n_next - n_ref;
        pot_half_next = exp(1i*k0*dn_eff_next*dz/2);

        A3 = A2 .* pot_half_next;

        % Apply boundary absorber
        E_rz(:,iz+1) = A3 .* absorber;
    end

    I_rz = abs(E_rz).^2;
end
