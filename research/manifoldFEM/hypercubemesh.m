function [node, isBdNode, isUpperBoundary] = hypercubemesh(dim, n)
%% HYPERCUBEMESH Generate node information for hypercube mesh
%
%   [node, isBdNode, isUpperBoundary] = hypercubemesh(dim, n) generates a 
%   structured mesh on a dim-dimensional hypercube with n element edges in each direction.
%
% Inputs:
%   dim - dimension of the hypercube
%   n   - number of elements (intervals) in each direction
%
% Outputs:
%   node           - Node coordinates matrix [n^dim x dim] with grid coordinates [0 to n] of the 1st node in each element
%   isBdNode       - Logical array [(n+1)^dim x 1], true for boundary nodes
%   isUpperBoundary - Logical array [(n+1)^dim x 1], true for upper boundary nodes
%   The upper boundary flags are used for stereographic projection for S^4
%   See also cubehexmesh, spherePoissonQ1
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

    numNode = (n+1)^dim;
    
    % Generate node coordinates using ndgrid (following cubehexmesh pattern)
    coords = cell(1, dim);
    grids = cell(1, dim);
    
    % Create coordinate arrays for each dimension
    coord_arrays = cell(1, dim);
    for i = 1:dim
        coord_arrays{i} = 0:n;
    end
    
    % Generate grid coordinates
    [grids{:}] = ndgrid(coord_arrays{:});
    for i = 1:dim
        coords{i} = grids{i}(:);
    end
    node = [coords{:}];
    
    % Initialize logical arrays for boundary flags
    isBdNode = false(numNode, 1);
    isUpperBoundary = false(numNode, 1);
    
    % Mark boundary and upper boundary nodes 
    for d = 1:dim
        % Nodes on lower boundary (coordinate = 0) or upper boundary (coordinate = N)
        bdNodes = (node(:, d) == 0) | (node(:, d) == n);
        isBdNode = isBdNode | bdNodes;
        
        % Nodes on upper boundary only (coordinate = n)
        upperNodes = (node(:, d) == n);
        isUpperBoundary = isUpperBoundary | upperNodes;
    end
end