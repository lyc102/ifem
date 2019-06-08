function [HB, NL, level] = HBstructure3(elem,N0,HBmesh)
%% HBSTRUCTURE3 reconstructs a hierachical structure of a 3-D mesh.
%
% HB(:,2:3) are two parent nodes of the node HB(:,1):
%             HB(:,2) --- HB(:,1) --- HB(:,3)
% NL records the range of indices in each level. The indices of nodes in
% the k-th level is given by NL(k)+1:NL(k+1).
%
% See also HBstructure3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('HBmesh','var'), HBmesh = []; end
N = max(elem(:));
HB = zeros(N,3);
% level = max(min(round(log2(N/N0)),16),2); % at least two level
level = 20;
NL(level+1) = N; 
for k = level: -1 : 2
    [elem,HBmesh,newHB] = uniformcoarsen3red(elem,HBmesh);  % try coasen red refinement
    if isempty(newHB) && ~isempty(HBmesh)     % then try coarsen bisection
        [elem,HBmesh,newHB] = uniformcoarsen3(elem,HBmesh);
    end
    if (isempty(newHB)) || (size(elem,1)< 2*N0) 
    % no nodes are removed or it reaches the coarsest level
        NL = NL(k:end);       
        break; 
    end
    NL(k) = NL(k+1) - size(newHB,1);
    HB(NL(k)+1:NL(k+1),1:3) = newHB(:,1:3);
end
level = length(NL)-1;