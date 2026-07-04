function n_rz = normalizeIndexMap_rz(n_complex, Nr, Nz)
%NORMALIZEINDEXMAP_RZ Convert scalar/vector/matrix n to size Nr x Nz.
%
% Accepted inputs:
%   scalar      -> constant medium
%   Nr x 1      -> radial profile replicated over z
%   1 x Nz      -> axial profile replicated over r
%   Nr x Nz     -> full map
%   Nz x Nr     -> transposed full map, automatically transposed

    if isscalar(n_complex)
        n_rz = repmat(n_complex, Nr, Nz);
        return;
    end

    s = size(n_complex);

    if isequal(s, [Nr, Nz])
        n_rz = n_complex;
    elseif isequal(s, [Nz, Nr])
        n_rz = n_complex.';
        warning('n_complex was Nz x Nr; transposed it to Nr x Nz.');
    elseif isvector(n_complex) && numel(n_complex) == Nr
        n_rz = repmat(n_complex(:), 1, Nz);
    elseif isvector(n_complex) && numel(n_complex) == Nz
        n_rz = repmat(reshape(n_complex, 1, Nz), Nr, 1);
    else
        error('n_complex must be scalar, Nr-vector, Nz-vector, Nr x Nz, or Nz x Nr.');
    end

    n_rz = double(n_rz);
end
