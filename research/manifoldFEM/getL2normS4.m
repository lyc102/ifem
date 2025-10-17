function L2 = getL2normS4(dim, elem, center, h, n, u)
%L2NORM Compute L2 norm using assembled mass matrix
%
%   L2 = L2Norm(dim, h, n, elem, vec)
%   computes the L2 norm of a vector on the d-dimensional sphere using 
%   the assembled mass matrix with volume form scaling.
%
% Inputs:
%   dim   - Dimension of the sphere (e.g., 4 for S^4)
%   h     - Mesh size in stereographic coordinates
%   n     - Number of elements in each direction
%   elem  - [numElem x 1] Global indices of first node of each element
%   u   - [(n+1)^dim x 1] FE solution compute norm of
%
% Outputs:
%   L2    - L2 norm: sqrt(vec' * M * vec)
%
% See also: assembleStiffness, getVolS4

    numNode = (n+1)^dim;
    numElem = length(elem);
    localNode = 2^dim;
    
    % Get local node mapping (same as in assembleStiffness)
    [localCoords, localNodeIdx] = getLocalNodeMapNdgrid(dim, n);
    
    % Compute local mass matrix integrals
    Integral = zeros(localNode, localNode);
    for i = 1:localNode
        for j = 1:localNode
            I = 1;
            for k = 1:dim
                if localCoords(i,k) == localCoords(j,k)
                    I = I * h/3;
                else
                    I = I * h/6;
                end
            end
            Integral(i,j) = I;
        end
    end
    vol = getVolS4(dim, center);
    
    % Assemble mass matrix
    nEntries = localNode * localNode;
    ii = zeros(nEntries*numElem, 1);
    jj = zeros(nEntries*numElem, 1);
    sM = zeros(nEntries*numElem, 1);
    
    index = 0;
    for i = 1:localNode
        for j = 1:localNode
            ii(index+1:index+numElem) = elem + localNodeIdx(i);
            jj(index+1:index+numElem) = elem + localNodeIdx(j);
            
            Mij = vol .* Integral(i,j);
            sM(index+1:index+numElem) = Mij;
            
            index = index + numElem;
        end
    end
    
    M = sparse(ii, jj, sM, numNode, numNode);
    
    % Compute L2 norm
    L2 = sqrt(u' * M * u);
end