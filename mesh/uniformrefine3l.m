function [node,elem,bdFlag,HB] = uniformrefine3l(node,elem,bdFlag)
%% UNIFORMREFINE3L uniformly refine a 3-D triangulation.
% 
% [node,elem] = uniformrefine3l(node,elem) divides each tetrahedra into
% eight small similar sub-tetrahedrons.
%
% [node,elem,bdFlag] = uniformrefine3l(node,elem,bdFlag) also update boundary
% conditions represented by bdFlag.
%
% [node,elem,~,HB] = uniformrefine3l(node,elem) outpus HB array which is
% useful for nodal interpolation.
%
% Difference with uniformrefine3: use the shortest diagonal.
%
% Example
%
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     elem = label3(node,elem);
%     figure(1); subplot(1,3,1); 
%     set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.5,0.3]);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%     [node,elem] = uniformrefine3(node,elem);
%     figure(1); subplot(1,3,2);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%     bdFlag = setboundary3(node,elem,'Dirichlet');
%     [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
%     figure(1); subplot(1,3,3);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%
% See also uniformbisect, uniformrefine, bisect, bisect3
% 
% Reference: S. Zhang. Successive subdivisions of tetrahedra and multigrid
% methods on t etrahedral meshes. Houston J. Math. 21, 541?556, 1995.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Construct data structure
[elem2dof,edge] = dof3P2(elem);
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Add new nodes
node(N+1:N+NE,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2; 
HB(:,[1 2 3]) = [(N+1:N+NE)', edge(:,1:2)]; 

%% Refine each tetrahedron into 8 tetrahedrons
t = 1:NT;
p(t,1:10) = elem2dof;
elem(8*NT,:) = [0 0 0 0]; % enlarge the elem array
elem(t,:)           = [p(t,1), p(t,5), p(t,6), p(t,7)];
elem(NT+1:2*NT,:)   = [p(t,5), p(t,2), p(t,8), p(t,9)];
elem(2*NT+1:3*NT,:) = [p(t,6), p(t,8), p(t,3), p(t,10)];
elem(3*NT+1:4*NT,:) = [p(t,7), p(t,9), p(t,10), p(t,4)];
% After cutting off four subtetrahedra at the corners, the remaining
% octahedron can be subdivided in three different ways by choosing one of
% three possible interior diagonals: [6 9], [7 8] or [5 10]. For better
% mesh quality, we choose the shortest interior diagonal. 
d1 = (node(p(t,7),1) - node(p(t,8),1)).^2 + ...
     (node(p(t,7),2) - node(p(t,8),2)).^2 + ...
     (node(p(t,7),3) - node(p(t,8),3)).^2;
d2 = (node(p(t,6),1) - node(p(t,9),1)).^2 + ...
     (node(p(t,6),2) - node(p(t,9),2)).^2 + ...
     (node(p(t,6),3) - node(p(t,9),3)).^2;
d3 = (node(p(t,5),1) - node(p(t,10),1)).^2 + ...
     (node(p(t,5),2) - node(p(t,10),2)).^2 + ...
     (node(p(t,5),3) - node(p(t,10),3)).^2;
[tempvar,minindex] = min([d1 d2 d3],[],2); %#ok<ASGLU>
% diagonal 7-8 is shorter than diagonal 6-9: use 7-8
idx = (minindex == 1);
elem(4*NT+t(idx),:) = [p(t(idx),6), p(t(idx),8), p(t(idx),7), p(t(idx),5)];
elem(5*NT+t(idx),:) = [p(t(idx),9), p(t(idx),7), p(t(idx),8), p(t(idx),5)];
elem(6*NT+t(idx),:) = [p(t(idx),6), p(t(idx),8), p(t(idx),10), p(t(idx),7)];
elem(7*NT+t(idx),:) = [p(t(idx),9), p(t(idx),7), p(t(idx),10), p(t(idx),8)];
% use diagonal 6-9
idx = (minindex == 2);
elem(4*NT+t(idx),:) = [p(t(idx),5), p(t(idx),6), p(t(idx),7), p(t(idx),9)];
elem(5*NT+t(idx),:) = [p(t(idx),5), p(t(idx),6), p(t(idx),9), p(t(idx),8)];
elem(6*NT+t(idx),:) = [p(t(idx),6), p(t(idx),7), p(t(idx),9), p(t(idx),10)];
elem(7*NT+t(idx),:) = [p(t(idx),6), p(t(idx),8), p(t(idx),10), p(t(idx),9)];
% use diagonal 5-10
idx = (minindex == 3);
elem(4*NT+t(idx),:) = [p(t(idx),6), p(t(idx),8), p(t(idx),10), p(t(idx),5)];
elem(5*NT+t(idx),:) = [p(t(idx),9), p(t(idx),7), p(t(idx),10), p(t(idx),5)];
elem(6*NT+t(idx),:) = [p(t(idx),6), p(t(idx),10), p(t(idx),5), p(t(idx),7)];
elem(7*NT+t(idx),:) = [p(t(idx),9), p(t(idx),10), p(t(idx),5), p(t(idx),8)];
% elem = label3(node,elem);

%% Update boundary edges
if (nargin>=3) && (~isempty(bdFlag))
    bdFlag(8*NT,:) = [0 0 0 0]; % enlarge the bdFlag array
    bdFlag(NT+1:2*NT,[1 3 4])   = bdFlag(t,[1 3 4]); 
    bdFlag(2*NT+1:3*NT,[1 2 4]) = bdFlag(t,[1 2 4]); 
    bdFlag(3*NT+1:4*NT,[1 2 3]) = bdFlag(t,[1 2 3]);
    idx = (minindex == 1);    % diagonal 7-8
    bdFlag(4*NT+t(idx),3) = bdFlag(t(idx),4);
    bdFlag(5*NT+t(idx),3) = bdFlag(t(idx),3);
    bdFlag(6*NT+t(idx),2) = bdFlag(t(idx),2);
    bdFlag(7*NT+t(idx),2) = bdFlag(t(idx),1);
    idx = (minindex == 2);   % diagonal 5-10
    bdFlag(4*NT+t(idx),2) = bdFlag(t(idx),3);
    bdFlag(5*NT+t(idx),3) = bdFlag(t(idx),4);
    bdFlag(6*NT+t(idx),3) = bdFlag(t(idx),2);
    bdFlag(7*NT+t(idx),1) = bdFlag(t(idx),1);
    idx = (minindex == 3);
    bdFlag(4*NT+t(idx),3) = bdFlag(t(idx),4);
    bdFlag(5*NT+t(idx),3) = bdFlag(t(idx),3);
    bdFlag(6*NT+t(idx),3) = bdFlag(t(idx),2);
    bdFlag(7*NT+t(idx),3) = bdFlag(t(idx),1);
    % change t in the last
    bdFlag(t,1) = 0;
else
    bdFlag = [];
end
