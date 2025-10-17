function idx = grid2ind(coords, base)
%% GRID2IND Convert grid coordinates to linear index in arbitrary base
%
%   idx = grid2ind(coords, base) converts multi-dimensional grid coordinates
%   to linear indices in a specified base system.
%
% Inputs:
%   coords - Grid coordinates [nPoints x dim], each row contains
%            coordinates (i1, i2, ..., i_dim) where ij ∈ [0, base-1]
%   base   - Base of the number system (grid size + 1 in each dimension)
%
% Outputs:
%   idx    - Linear index (0-indexed) [nPoints x 1]
%
% Example:
%   % Convert 4D grid coordinates to node index in 16x16x16x16 mesh
%   idx = grid2ind([0,0,0,0], 17);    % Returns 0
%   idx = grid2ind([1,2,3,4], 17);    % Returns linear index
%
%   % Verify round-trip conversion
%   coords = [1, 2, 3, 4];
%   idx = grid2ind(coords, 17);
%   coords_back = ind2grid(idx, 17, 4);  % Should equal original coords
%
% See also: ind2grid, sub2ind, ind2sub

    dim = size(coords, 2);
    idx = 0;
    
    for i = 1:dim
        idx = idx + coords(:, i) * base^(dim - i);
    end
end


  % Columns 1 to dim: Local grid coordinates within a hypercube element

    % Each entry is either 0 or 1
    % Represents the relative position of a local node within the reference element
    % For example, in 4D:
    