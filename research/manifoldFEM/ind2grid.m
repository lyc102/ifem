function coords = ind2grid(idx, base, dim)
%% IND2GRID Convert linear index to grid coordinates in arbitrary base
%
%   coords = ind2grid(idx, base, dim) converts linear indices to 
%   multi-dimensional grid coordinates in a specified base system.
%
% Inputs:
%   idx   - Linear index (or indices) to convert (0-indexed)
%   base  - Base of the number system (grid size + 1 in each dimension)
%   dim   - Number of dimensions
%
% Outputs:
%   coords - Grid coordinates [numel(idx) x dim], each row contains
%            coordinates (i1, i2, ..., i_dim) where ij ∈ [0, base-1]
%
% Example:
%   % Convert node index to 4D grid coordinates in 16x16x16x16 mesh
%   coords = ind2grid(0, 17, 4);    % Returns [0,0,0,0]
%   coords = ind2grid(100, 17, 4);  % Returns grid position of node 100
%
% See also: ind2sub, sub2ind

    nIndices = numel(idx);
    coords = zeros(nIndices, dim);
    
    for i = 1:dim
        coords(:, i) = floor(idx / base^(dim - i));
        idx = idx - coords(:, i) * base^(dim - i);
    end
end