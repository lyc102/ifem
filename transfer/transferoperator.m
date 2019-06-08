function [Pro,Res,freeNodej] = transferoperator(HB,NL,isFreeNode)
%% TRANSFEROPERATOR transfer operators between multilevel meshes
%
%  [Pro,Res] = transferoperator(HB,NL) construct matrix representation of
%  transfer operators between levels. All output are cell(level,1).
%  
%  In the input, HB is the hierarchical structure of nodes. NL records the
%  starting and ending indices of nodes in each level. see HBstructure or
%  HBstructure3.
% 
%  Prolongation (Pro) and restrction (Res) is constructed using standard
%  linear finite element basis and Res = Pro'. 
%
%  See also mg, HBstructure, HBstructure3.
%
%  Reference: L. Chen. MultiGrid Methods.

if ~exist('isFreeNode','var'), isFreeNode = []; end 
level = length(NL)-1; 
Pro = cell(level,1);
Res = cell(level,1);
freeNodej = cell(level,1);
for j = level:-1:2
    % free nodes in each level
    freeNodej{j} = find(isFreeNode);
    % fine node and coarse node index
    fineNodeRange = NL(j)+1:NL(j+1);
    fineNode = HB(fineNodeRange,1);
    nFineNode = NL(j+1)-NL(j);
    coarseNode = (1:NL(j))';
    isCoarseNode = true(NL(j+1),1);
    isCoarseNode(fineNode) = false;
    coarseNodeFineIdx = find(isCoarseNode);
    ii = [coarseNodeFineIdx; fineNode; fineNode];
    jj = [coarseNode; HB(fineNodeRange,2); HB(fineNodeRange,3)];
    ss = [ones(NL(j),1); 0.5*ones(nFineNode,1); 0.5*ones(nFineNode,1)];
    % remove fixed node from the index
    if ~isempty(isFreeNode)
        isFreeNodeC = isFreeNode;
        isFreeNodeC(fineNode) = [];    
        Nfree = sum(isFreeNode);
        NfreeC = sum(isFreeNodeC);
        idx = isFreeNode(ii) & isFreeNodeC(jj);
        idxMap(isFreeNode) = (1:Nfree)';
        idxMapC(isFreeNodeC) = (1:NfreeC)';
        ii = idxMap(ii(idx));
        jj = idxMapC(jj(idx));
        ss = ss(idx);
        % generate prolongation matrix
        Pro{j-1} = sparse(ii,jj,ss,Nfree,NfreeC);
        isFreeNode = isFreeNodeC;    
    else
        Pro{j-1} = sparse(ii,jj,ss,NL(j+1),NL(j));
    end
    Res{j} = Pro{j-1}';                       
end