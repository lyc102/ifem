function rhs = assembleRHS(pde, center, elem, h, n, vol, chart)
%ASSEMBLERHS Assemble right-hand side vector for Laplace-Beltrami equation on S^d
%
%   rhs = assembleRHS(h, n, center, elem, vol, pde, chart)
%   assembles the global right-hand side vector for the Laplace-Beltrami equation
%   on the d-dimensional sphere using Q1 finite elements in stereographic 
%   projection coordinates.
%
%   The right-hand side corresponds to:
%       ∫_{S^d} f * v * dV_g
%
%   where f is the source term, v is the test function, and dV_g is the 
%   volume form on S^d.
%
% Inputs:
%   pde.dim              - Embedding dimension (sphere is S^dim)
%   pde.b                - L2 coefficient (scalar)
%   h                - Mesh size in stereographic coordinates
%   n                - Number of elements in each direction
%   elem - [numElem x 1] Global indices of first node of each element
%                      (elements that don't touch upper boundary)
%   center           - [numElem x dim] Element centers in stereographic coordinates
%   vol              - [numElem x 1] Volume form values at element centers
%   pde              - PDE data structure with fields:
%                      .f: function handle f(center, b, chart) returning source term
%   chart            - 'upper' or 'lower' stereographic projection chart
%
% Outputs:
%   rhs - [(n+1)^dim x 1] Right-hand side vector
%         rhs(i) = ∫ f * \phi_i * vol
%
% Algorithm:
%   1. Evaluate source term f at element centers (piecewise constant approximation)
%   2. Compute local contributions using tensor product quadrature
%   3. Scale by volume form at element centers
%   4. Assemble global vector using local-to-global index mapping (vectorized)
%
% Note:
%   - Uses Q1 (multilinear) basis functions on hypercube elements
%   - Source term f is approximated as piecewise constant on each element
%   - Each element has 2^dim local nodes
%   - Volume form: 2^dim * (1 + |x|^2)^(-dim)
%
% See also: assembleStiffness, getVolS4, getLocalNodeMapNdgrid
    b = pde.b;
    dim = pde.dim;
    numNode = (n+1)^dim;
    numElem = length(elem);
    localNode = 2^dim;
    
    % Local node relative positions using ndgrid
    [~, localNodeIdx] = getLocalNodeMapNdgrid(dim, n);
    
    % Get f function values at element centers (piecewise constant)
    f = pde.f(center, b, chart);
    
    % Each local basis function integrates to (h/2)^dim over the element
    localVol = (h * 0.5)^dim;
    
    % indices and values for all element-local node pairs
    ii = zeros(numElem * localNode, 1);
    vv = zeros(numElem * localNode, 1);
    
    for s = 1:localNode
        idx = (s-1)*numElem + 1 : s*numElem;
        
        % Global node indices for all elements
        ii(idx) = elem + localNodeIdx(s);
        vv(idx) = vol .* f * localVol;
    end
    
    rhs = accumarray(ii, vv, [numNode, 1]);
end