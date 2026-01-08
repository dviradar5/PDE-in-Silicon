function L = lap1dNeumann(N, d)
% L*u ≈ d^2u/dx^2 with Neumann BC (du/dx=0) at both ends.

    e = ones(N,1);
    S = spdiags([e -2*e e], [-1 0 1], N, N);

    % Neumann BC via ghost reflection:
    S(1,2)   = 2;
    S(N,N-1) = 2;

    L = S / d^2;
end
