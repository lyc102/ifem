function Pro = transferelem(elemc,elemf,tree)
%% TRANSFERELEM constructs prolongation operator 
%
% Pro = TRANSFERELEM(elemc,elemf,tree) constructs the prolongation matrix
% Pro from the coarse mesh elemc to the fine mesh elemf. The tree array
% records the parent tree(:,1) and its two children tree(:,2:3).
%
% Example
%
%   [node,elem] = squaremesh([0 1 0 1],1/2);
%   p = (1:size(elem,1))'; 
%   figure(1);
%   subplot(1,3,1); showsolution(node,elem,p,'EdgeColor','k'); view(2);
%   [nodef,elemf,~,~,tree] = bisect(node,elem,[1 2]);
%   Pro = transferelem(elem,elemf,tree);
%   newp = Pro*p;
%   subplot(1,3,2); showsolution(nodef,elemf,newp,'EdgeColor','k'); view(2);
%   p2 = eleminterpolate(p,tree);
%   display(norm(newp - p2));
%   node = nodef; elem = elemf;
%   [nodef,elemf,~,~,tree] = bisect(nodef,elemf,[1 2]);
%   Pro = transferelem(elem,elemf,tree);
%   newp = Pro*p2;
%   p3 = eleminterpolate(p2,tree);
%   subplot(1,3,3); showsolution(nodef,elemf,newp,'EdgeColor','k'); view(2);
%   display(norm(newp - p3));
%
%
% See also eleminterpolate, transferedge

% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%%
NTc = size(elemc,1);
NTf = size(elemf,1);
inCoarse = true(NTf,1);
inCoarse(tree(:,3)) = false; % always remove R
% coarseElemFineIdx = find(inCoarse); 
ii = (1:NTf)';
jj = zeros(NTf,1);
% for each element in the fine mesh, find one coarse element
% Case 0: coarse element
jj(inCoarse) = (1:NTc)';
% Case 1: The parent element is in the coarse mesh.
idx1 = (tree(:,1) <= NTc);
jj(tree(idx1,2)) = tree(idx1,1);
jj(tree(idx1,3)) = tree(idx1,1);
% Case 2: The parent is in the intermediate mesh. A triangle t could
% be bisected twice and a child could be a parent
parent(tree(idx1,3)) = tree(idx1,1);
idx2 = (tree(:,1) > NTc);
grandParent = parent(tree(idx2,1));
jj(tree(idx2,2)) = grandParent;
jj(tree(idx2,3)) = grandParent;
Pro = sparse(ii,jj,1,NTf,NTc);