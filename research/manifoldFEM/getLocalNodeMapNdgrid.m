function [localCoords, localNodeIdx] = getLocalNodeMapNdgrid(dim, n)
%GETLOCALNODEMAPNDGRID Generate local to global node index mapping for hypercube elements
%
%   [localCoords, localNodeIdx] = getLocalNodeMapNdgrid(dim, n)
%
% Inputs:
%   dim - Dimension of the hypercube
%   n   - Number of elements in each direction
%
% Outputs:
%   localCoords - [2^dim x dim] Local grid binary-valued coordinates in {0,1}^dim 
%   Each row represents a corner of the reference hypercube
%   localNodeIdx    - [2^dim x 1] Global node index offset
%                 Offset to add to element's first node to get global index
%  Example in 4D:
%  [0,0,0,0] = first corner (origin of the element)
%  [1,0,0,0] = corner in the +x_1 direction
%  [1,1,1,1] = opposite corner (furthest from origin)

%  The offset to add to the element's first node index (nodeIdxNonOverlap) to get this local node's global (linear) index
%  Computed using the base-(n+1) number system
%  For example, if element starts at global node 100:
%  Local node with offset 0 → global node 100
%  Local node with offset 5 → global node 105

    % Generate all combinations of {0,1}^dim using ndgrid
    coords = cell(1, dim);
    grids = cell(1, dim);
    [grids{:}] = ndgrid([0, 1]);
    
    for i = 1:dim
        coords{i} = grids{i}(:);
    end
    
    localCoords = [coords{:}];
    
    % Compute global index offset for each local node
    localNodeIdx = grid2ind(localCoords, n+1);
end