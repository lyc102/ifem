function u = nodeinterpolate(u,HB)
%% NODEINTERPOLATE interpolate a piecewise linear function.
%
% u = nodeinterpolate(u,HB) interpolate a linear finite element function u
% from a coarse grid to a fine grid or the other way around. The input
% array HB(:,1:3) records hierarchical structure of new nodes HB(:,1) added
% from a coarse grid to a fine one. HB(:,2:3) are two parent nodes of
% HB(:,1). It can be obtained by bisection. Going from fine to coarse, HB =
% indexMap which is the map between the indices of nodes in the fine grid
% to that of the coarse grid. For example, indexMap(10) = 6 means the 10-th
% node in the fine grid is now the 6-th node in the coarse one. Therefore
% indexMap(k) = 0 means k is removed. indexMap is obtained by coarsening.
%
% Example 1: 2D bisect and coarsen
%   node = [0,0; 1,0; 1,1; 0,1];
%   elem = [2,3,1; 4,1,3];      
%   u = [0 0 0 1]; 
%   figure(1);
%   subplot(1,3,1); showsolution(node,elem,u,'EdgeColor','k');
%   [node,elem,~,HB] = bisect(node,elem,1);
%   u = nodeinterpolate(u,HB);
%   subplot(1,3,2); showsolution(node,elem,u,'EdgeColor','k');
%   [node,elem,~,~,indexMap] = coarsen(node,elem,1:size(elem,1));
%   u = nodeinterpolate(u,indexMap);
%   subplot(1,3,3); showsolution(node,elem,u,'EdgeColor','k');
%
% Example 2: 3D bisect and coarsen
%     [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],2);
%     u = ones(size(node,1),1);
%     [node,elem,~,HB] = bisect3(node,elem,[1 2],[],HB);
%     [node,elem,~,HB] = bisect3(node,elem,[1 2],[],HB);
%     [node,elem,~,HB] = bisect3(node,elem,[1 2],[],HB);
%     showmesh3(node,elem);
%     u = nodeinterpolate(u,HB);
%     [node,elem,~,HB,indexMap] = coarsen3(node,elem,'all',[],HB);
%     u = nodeinterpolate(u,indexMap);
%
% See also bisect, coarsen, eleminterpolate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%%
if isempty(HB), return; end
if size(u,2) > size(u,1) % u is a column vector
    u = u';
end
oldN = size(u,1);
newN = max(size(HB,1),max(HB(:,1)));
if oldN >= newN % fine grid to coarse grid
    idx = (HB == 0);  % HB is indexMap output by coarsening
    u(idx,:) = [];
else            % coarse grid to fine grid
    if min(HB(:,1))>oldN % only new nodes are recorded in HB (2-D bisection)
        u(HB(1:end,1),:) = (u(HB(1:end,2),:)+u(HB(1:end,3),:))/2;
    else        % new nodes is stored starting from oldN+1 (3-D bisection)
        u(newN) = 0;  % preallocation
        for k = oldN+1:newN  % have to perform the interpolation sequentially
            u(HB(k,1),:) = (u(HB(k,2),:) + u(HB(k,3),:))/2;            
        end
    end
end