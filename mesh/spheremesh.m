function [node,elem] = spheremesh(N)
%% SPHEREMESH generate uniform mesh on unit sphere surface.
%
% Author: Huayi Wei <huayiwei1984@gmail.com> 

surfacedata = spheresurface();
[node,elem] = surfacedata.initmesh();

for i = 1:N
   [node,elem] = smeshuniformrefine(node,elem);
   node = surfacedata.project(node);
   optsurfacemesh(node,elem,surfacedata);
end