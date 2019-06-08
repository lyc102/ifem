function newp = eleminterpolate(p,tree,weight)
%% ELEMINTERPOLATE interpolate a piecewise constant function. 
%
% newp = eleminterpolate(p,tree) interpolate a piecewise constant function
% p from a coarse grid to a fine grid or the other way around. tree(:,1:3)
% stores the binary tree of the coarsening. tree(:,1) is the index of
% parent element in the coarsened mesh and tree(:,2:3) are two children
% indices in the original mesh.
%
% For uniform refinement, the prolongation from the coarse grid to the fine
% grid is simply |newp = repmat(p,4,1)| and from the fine grid to the
% coarse one is |newp = p(1:NT/4)|.
%
% Example
%
%   [node,elem] = squaremesh([0 1 0 1],1/2);
%   p = 1:size(elem,1); 
%   figure(1);
%   subplot(1,4,1); showsolution(node,elem,p,'EdgeColor','k'); view(2);
%   [node,elem,~,~,tree] = bisect(node,elem,[1 2]);
%   p = eleminterpolate(p,tree);
%   subplot(1,4,2); showsolution(node,elem,p,'EdgeColor','k'); view(2);
%   [node,elem,~,~,tree] = bisect(node,elem,[1 2]);
%   p = eleminterpolate(p,tree);
%   subplot(1,4,3); showsolution(node,elem,p,'EdgeColor','k'); view(2);
%   [node,elem,~,~,tree] = coarsen(node,elem,'all');
%   p = eleminterpolate(p,tree);
%   subplot(1,4,4); showsolution(node,elem,p,'EdgeColor','k'); view(2);
%   
% See also bisect, coarsen, nodeinterpolate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if ~exist('weight','var'), weight = 1; end
%%
newp = p;
if (size(tree,1)==0), return; end % no change
NTin = length(p); 
NTf = max(tree(:,3)); 
if NTin < NTf 		   % coarse grid to fine grid    
    NTc = NTin;
    newp = zeros(NTf,1);
    isNew = tree(:,3); % Right child is a new element
    inCoarse = true(NTf,1);
    inCoarse(isNew) = false;
    % Case 0: still in coarse element
    newp(inCoarse) = p;
    % Case 1: The parent element is in the coarse mesh. A triangle t could
    % be bisected twice and a child could be a parent
    idx = (tree(:,1) <= NTc);
    if weight == 1
        newp(tree(idx,2)) = p(tree(idx,1));
        newp(tree(idx,3)) = p(tree(idx,1));
    else
        newp(tree(idx,2)) = weight*p(tree(idx,1));
        newp(tree(idx,3)) = weight*p(tree(idx,1));        
    end
    % Case 2: The parent is in the intermediate mesh
    idx = (tree(:,1) > NTc);
    if weight == 1
        newp(tree(idx,2)) = newp(tree(idx,1));
        newp(tree(idx,3)) = newp(tree(idx,1));
    else
        newp(tree(idx,2)) = weight*newp(tree(idx,1));
        newp(tree(idx,3)) = weight*newp(tree(idx,1));
    end        
elseif NTin == NTf      % fine grid to coarse grid
    newp = p;
    newp(tree(:,3)) = [];
end