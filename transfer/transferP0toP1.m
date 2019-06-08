function [Pro,Res] = transferP0toP1(node,elem,fixedNode)
%% TRANSFERP0TOP1: the integral average transfer from P0 to P1
%
%  Created by Ming Wang, at March, 2012. 
%
N = size(node,1);
NT = size(elem,1);
[~,area] = gradbasis(node,elem);
patchArea = accumarray(elem(:),repmat(area,3,1),[N 1]); 
invPatchArea = spdiags(1./patchArea,0,N,N);
pro = sparse(elem(:),repmat((1:NT)',3,1),repmat(area,3,1),N,NT);
Pro = invPatchArea*pro;
if exist('fixedNode','var'), Pro(fixedNode,:) = 0; end
Res = Pro';