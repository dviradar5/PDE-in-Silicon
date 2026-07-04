function mask = radialAbsorbingMask(r, a, b)
%RADIALABSORBINGMASK Smooth multiplicative sponge at the radial edge.
%
% mask = exp(-a*(r/rmax)^b). In your previous BPM code this kind of mask is
% often used only to suppress numerical edge reflections/sidelobes. Set a=0
% or [] to disable.

    if nargin < 2 || isempty(a), a = 0; end
    if nargin < 3 || isempty(b), b = 8; end

    r = r(:);
    if a == 0
        mask = ones(size(r));
        return;
    end

    rmax = max(abs(r));
    if rmax == 0
        mask = ones(size(r));
    else
        mask = exp(-a * (abs(r) / rmax).^b);
    end
end
