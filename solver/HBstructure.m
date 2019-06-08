function [HB, NL, level] = HBstructure(elem,N0)
%% HBSTRUCTURE reconstruct a hierachical structure of the mesh.
%
% HB(:,2:3) are two parent nodes of the node HB(:,1):
%             HB(:,2) --- HB(:,1) --- HB(:,3)
% NL records the range of indices in each level. The indices of nodes in
% the k-th level is given by NL(k)+1:NL(k+1).
%
% See also HBstructure3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = max(elem(:));
HB = zeros(N,3);
% level = max(min(round(log2(N/N0)),16),2); % at least two level
level = 20;
NL(level+1) = N; % now NL(1:level) = 0;
for k = level: -1 : 2
    [elem,newHB] = uniformcoarsenred(elem);  % try coasen red refinement first
    if isempty(newHB) && isequal(size(elem,2),3)
        [elem,newHB] = uniformcoarsen(elem); % then coarsen bisection
    end
    if isempty(newHB) && isequal(size(elem,2),4)
        [elem,newHB] = uniformcoarsenquadred(elem); % coarsen uniform quad mesh
    end
    if ~isempty(newHB) 
        NL(k) = NL(k+1) - size(newHB,1); % update NL(k)
        HB(NL(k)+1:NL(k+1),1:3) = newHB(:,1:3);
    else   % no nodes are removed
        NL = NL(k:end);       
        break; 
    end
    if  size(elem,1)< 2*N0 % it reaches the coarsest level
        NL = NL(k-1:end);
        break
    end
end
level = length(NL)-1;    % actual level