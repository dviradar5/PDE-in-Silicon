function L = lap1dNeumannCylR(r, dr)

    r = r(:);
    Nr = numel(r);

    main  = (-2/dr^2) * ones(Nr,1);
    upper = zeros(Nr-1,1);
    lower = zeros(Nr-1,1);

    % Interior nodes i = 2..Nr-1
    for i = 2:Nr-1
        ri = r(i);
        lower(i-1) = 1/dr^2 - 1/(2*dr*ri);
        upper(i)   = 1/dr^2 + 1/(2*dr*ri);   % note: upper is Nr-1, so index i<=Nr-1 is OK
    end

    % r = 0 (i=1): L p ≈ 4 (p2 - p1)/dr^2
    main(1)  = -4/dr^2;
    upper(1) =  4/dr^2;

    % r = R (i=Nr): Neumann dp/dr=0 -> L p ≈ 2 (p_{N-1} - p_N)/dr^2
    main(Nr)     = -2/dr^2;
    lower(Nr-1)  =  2/dr^2;

    % Pad to Nr so spdiags can concatenate
    lower_full = [lower; 0];
    upper_full = [0; upper];

    L = spdiags([lower_full main upper_full], [-1 0 1], Nr, Nr);
end