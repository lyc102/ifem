function [node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB)
%% UNIFORMBISECT3 uniformly bisect a 3-D triangulation. 
%
% [node,elem] = uniformbisect3(node,elem) divides each tetrahedron into 8 small
% tetrahedrons using longest vertex bisection.
%
% [node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB) update HB and
% boundary conditions.
%
% Example
%
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     elem = label3(node,elem);
%     figure(1); subplot(1,3,1); 
%     set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.5,0.3]);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%     [node,elem] = uniformbisect3(node,elem);
%     figure(1); subplot(1,3,2);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%     bdFlag = setboundary3(node,elem,'Dirichlet');
%     [node,elem,bdFlag] = uniformbisect3(node,elem,bdFlag);
%     figure(1); subplot(1,3,3);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%
% See also uniformbisect, uniformrefine, bisect, bisect3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('bdFlag','var')
    bdFlag =[]; 
end
if ~exist('HB','var')
    HB = []; 
end
[node,elem,bdFlag,HB] = bisect3(node,elem,'all',bdFlag,HB);
[node,elem,bdFlag,HB] = bisect3(node,elem,'all',bdFlag,HB);
[node,elem,bdFlag,HB] = bisect3(node,elem,'all',bdFlag,HB);