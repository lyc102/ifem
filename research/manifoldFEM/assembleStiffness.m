function A = assembleStiffness(pde, center, elem, h, n, vol)
%ASSEMBLESTIFFNESS Assemble stiffness matrix for Laplace-Beltrami operator on S^d
%
%   A = assembleStiffness(pde, center, elem, h, n, vol)
%   assembles the subdomain stiffness matrix for the Laplace-Beltrami operator
%   with L2 mass term on the d-dimensional sphere using tensorial finite elements on each subdomain using stereographic projection coordinates.
%
%   The weak form is:
%       ∫_{S^d} \nabla_g u \cdot \nabla_g v + b*u*v dV
%
%   where \nabla_g is the "surface" gradient and dV is the volume form on S^d.
%
% Inputs:
%   pde.dim              - 4 (sphere is S^dim)
%   pde.b                - L2 coefficient (scalar)
%   elem - [numElem x 1] Global indices of first node of each element, in domain decomposition it is elemNonOverlap
%                      (elements that do not touch upper boundary)
%   h                - Mesh size in stereographic coordinates
%   n                - Number of elements in each direction
%   center           - [numElem x dim] Element centers in stereographic coordinates
%   vol              - [numElem x 1] Volume form values at element centers
%
% Outputs:
%   A - [(n+1)^dim x (n+1)^dim] Sparse stiffness matrix
%       A(i,j) = \int g^{kl} \partial_{x^k} \phi_i \partial_{x^l} \phi_j + b*\phi_i*\phi_j * vol
%       where g^{kl} is the inverse metric tensor from stereographic projection
%
% Algorithm:
%   1. Compute local element stiffness matrices using tensor product quadrature
%   2. Scale by metric tensor coefficient for Laplace-Beltrami operator
%   3. Add L2 mass term scaled by volume form
%   4. Assemble global matrix using local-to-global index mapping
%
% Note:
%   - Uses multilinear basis functions on hypercube elements
%   - Each element has 2^dim local nodes
%   - Metric coefficient: 2^(dim-2) * (1 + |x|^2)^(2-dim)
%   - Volume form: 2^dim * (1 + |x|^2)^(-dim)
%
% Example:
%   dim = 4; n = 8; b = 1;
%   [node, isBdNode, isUpperBoundary] = hypercubemesh(dim, n);
%   elemNonOverlap = find(~isUpperBoundary);
%   r = 1.15; h = 2*r/n;
%   nodeInt = node(~isUpperBoundary, :);
%   center = -r + (nodeInt + 0.5) * h;
%   vol = getVolS4(dim, center);
%   A = assembleStiffness(dim, elemNonOverlap, b, h, n, center, vol);
%
% See also: getMetricS4, getVolS4
%
% Reference:
%   Cao, Shuhao, and Lizhen Qin. "A numerical domain decomposition method 
%   for solving elliptic equations on manifolds." SIAM Journal on 
%   Scientific Computing 46, no. 1 (2024): A376-A398.

    dim = pde.dim;
    b = pde.b;
    numNode = (n+1)^dim;
    numElem = n^dim;
    localNode = 2^dim;
    
    % Local node relative positions using ndgrid
    [localCoords, localNodeIdx] = getLocalNodeMapNdgrid(dim, n);
    
    % Local integrals
    Integral = zeros(localNode, localNode);
    DIntegral = zeros(localNode, localNode);
    for i = 1:localNode
        for j = 1:localNode
            I = 1;
            DIoverI = 0;
            for k = 1:dim
                if localCoords(i,k) == localCoords(j,k)
                    I = I * h/3;
                    DIoverI = DIoverI + 3/(h^2);
                else
                    I = I * h/6;
                    DIoverI = DIoverI - 6/(h^2);
                end
            end
            Integral(i,j) = I;
            DIntegral(i,j) = I * DIoverI;
        end
    end
    
    % Second order term (metric-dependent coefficient for Laplace-Beltrami)
    secondOrderCoeff = getMetricS4(dim, center);
    
    nEntries = localNode * localNode;
    ii = zeros(nEntries*numElem, 1); 
    jj = zeros(nEntries*numElem, 1);
    sA = zeros(nEntries*numElem, 1);
    
    index = 0;
    for i = 1:localNode
        for j = 1:localNode
            % Local to global index map
            ii(index+1:index+numElem) = elem + localNodeIdx(i);
            jj(index+1:index+numElem) = elem + localNodeIdx(j);
            
            % Element stiffness entry (gradient term + mass term)
            Aij = secondOrderCoeff .* DIntegral(i,j) + b * vol .* Integral(i,j);
            sA(index+1:index+numElem) = Aij;
            
            index = index + numElem;
        end
    end
     
    A = sparse(ii, jj, sA, numNode, numNode);
end