function [bdNodeIdxOther, bdInterpWeights] = transitionChartMap(node, bdNode, r, n)
%TRANSITIONCHARTMAP Compute boundary interpolation data for domain decomposition
%
%   [bdNodeIdxOther, bdInterpWeights] = transitionChartMap(dim, node, bdNode, r, n)
%   computes the interpolation data needed to transfer boundary values between
%   overlapping stereographic charts in the Schwarz domain decomposition method.
%
% Inputs:
%   dim     - Embedding dimension (sphere is S^dim)
%   node    - [numNode x dim] Node coordinates in grid coordinates [0, n]
%   bdNode  - [numBdNode x 1] Indices of boundary nodes in current chart
%   r       - Overlapping radius in stereographic coordinates
%   n       - Number of elements in each direction
%
% Outputs:
%   bdNodeIdxOther  - [numBdNode x 2^dim] Global node indices in opposite chart
%                     Each row contains the 2^dim nodes of the element that
%                     contains the stereographically inverted boundary point
%   bdInterpWeights - [numBdNode x 2^dim] Interpolation weights for each boundary node
%                     Sum of each row equals 1
%
% Algorithm:
%   1. Map boundary nodes to opposite chart via stereographic inversion
%   2. Find which element contains each mapped point
%   3. Compute interpolation weights within that element
%
% Mathematical Details:
%   - Stereographic inversion: x_other \approx r^2 * x_current / |x_current|^2
%   - Continuous coordinates Theta = (r + x/|x|^2) / h map to grid [0,n]
%   - Integer part floor(Theta) identifies the containing element
%   - Fractional part (Theta - floor(Theta)) determines interpolation weights
%   - Multilinear weight: \prod_i (1 - |\theta_i - \delta_i|) where \delta_i ∈ {0,1}
%
% Example:
%   dim = 4; n = 8; r = 1.15;
%   [node, isBdNode, ~] = hypercubemesh(dim, n);
%   bdNode = find(isBdNode);
%   [bdIdx, bdWeights] = transitionChartMap(dim, node, bdNode, r, n);
%   % Transfer values: u_bd = bdWeights .* u_other(bdIdx)
%
% See also: getLocalNodeMapNdgrid

    dim = size(node, 2);
    h = 2*r/n;
    numBdNode = size(bdNode, 1);
    localNode = 2^dim;
    
    % Convert boundary nodes to stereographic coordinates
    pt = -r + node(bdNode, :) * h;
    
    % Map to opposite chart via stereographic inversion
    % Formula: x_other = r + x / |x|²
    normSq = dot(pt, pt, 2);
    ptOther = r + pt ./ normSq;
    
    % Convert to continuous grid coordinates in opposite chart
    Theta = ptOther / h;
    ThetaInt = floor(Theta);
    ThetaFrac = Theta - ThetaInt;  % Fractional part in [0,1]^dim
    
    % Find first node of containing element in opposite chart
    firstNodeIdx = grid2ind(ThetaInt, n+1) + 1;  % +1 for 1-based indexing
    
    % Get local node mapping
    [localCoords, localNodeOffset] = getLocalNodeMapNdgrid(dim, n);
    
    % Preallocate output arrays
    bdNodeIdxOther = zeros(numBdNode, localNode);
    bdInterpWeights = zeros(numBdNode, localNode);
    
    % Compute global indices and interpolation weights for each local node
    for j = 1:localNode
        % Global node index in opposite chart
        bdNodeIdxOther(:, j) = firstNodeIdx + localNodeOffset(j);
        
        weight = ones(numBdNode, 1);
        for s = 1:dim
            % Distance from local node corner in dimension s
            dist = abs(ThetaFrac(:, s) - localCoords(j, s));
            weight = weight .* (1 - dist);
        end
        bdInterpWeights(:, j) = weight;
    end
end