function [node,elem,bdFlag,HB] = uniformrefine3(node,elem,bdFlag,HB)
%% UNIFORMREFINE3 uniformly refine a 3-D triangulation.
% 
% [node,elem] = uniformrefine3(node,elem) divides each tetrahedra into
% eight small similar sub-tetrahedrons.
%
% [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag) also update boundary
% conditions represented by bdFlag.
%
% [node,elem,~,HB] = uniformrefine3(node,elem) outpus HB array which is
% useful for nodal interpolation. Unlike bisect3, HB is not useful for the
% coarsening. See uniformcoarsen3red.
%
% Warning: uniform refine in 3D is not orientation presereved. Some volumes
% of the sub-tetrahedron are negative. The ordering is special such that
% the refinement can be easy. See uniformrefine3doc for details.
%
% Example
%
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1 2 3 7; 1 4 3 7; 1 5 6 7; 1 5 8 7; 1 2 6 7; 1 4 8 7];
%     figure(1); subplot(1,3,1); 
%     set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.5,0.3]);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%     [node,elem] = uniformrefine3(node,elem);
%     figure(1); subplot(1,3,2);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%     bdFlag = setboundary3(node,elem,'Dirichlet');
%     [node,elem,~,bdFlag] = uniformrefine3(node,elem,[],bdFlag);
%     figure(1); subplot(1,3,3);
%     showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%
% See also uniformbisect, uniformrefine, bisect, bisect3, uniformcoarsen3red
%
% Reference page in Help browser
%  <a href="matlab:ifem uniformrefine3doc">ifem uniformrefine3doc</a>
%
% Reference: J. Bey. Simplicial grid refinement: on Freudenthal's algorithm
% and the optimal number of congruence classes. Numer. Math.. 85(1):1--29,
% 2000. p11 Algorithm: RedRefinement3D.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('bdFlag','var'), bdFlag =[]; end
if ~exist('HB','var'),    HB = [];     end

%% Construct data structure
[elem2dof,edge] = dof3P2(elem);
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Add new nodes
node(N+1:N+NE,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2;
if ~isempty(HB)    
    maxgeneration = max(HB(:,4));
    HB(N+1:N+NE,[1 2 3]) = [(N+1:N+NE)', edge(:,1:2)]; 
    HB(N+1:N+NE,4) = maxgeneration + 1; 
end

%% Refine each tetrahedron into 8 tetrahedrons
t = 1:NT;
p(t,1:10) = elem2dof;
elem(8*NT,:) = [0 0 0 0]; % enlarge the elem array
elem(t,:) = [p(t,1), p(t,5), p(t,6), p(t,7)];
elem(NT+1:2*NT,:) = [p(t,5), p(t,2), p(t,8), p(t,9)];
elem(2*NT+1:3*NT,:) = [p(t,6), p(t,8), p(t,3), p(t,10)];
elem(3*NT+1:4*NT,:) = [p(t,7), p(t,9), p(t,10), p(t,4)];
% always use diagonal 6-9. The ordering is important. See the reference.
elem(4*NT+1:5*NT,:) = [p(t,5), p(t,6), p(t,7), p(t,9)];
elem(5*NT+1:6*NT,:) = [p(t,5), p(t,6), p(t,8), p(t,9)];
elem(6*NT+1:7*NT,:) = [p(t,6), p(t,7), p(t,9), p(t,10)];
elem(7*NT+1:8*NT,:) = [p(t,6), p(t,8), p(t,9), p(t,10)];

%% Update boundary edges
if ~isempty(bdFlag)
    bdFlag(8*NT,:) = [0 0 0 0]; % enlarge the bdFlag array
    bdFlag(NT+1:2*NT,[1 3 4]) = bdFlag(t,[1 3 4]); 
    bdFlag(2*NT+1:3*NT,[1 2 4]) = bdFlag(t,[1 2 4]); 
    bdFlag(3*NT+1:4*NT,[1 2 3]) = bdFlag(t,[1 2 3]);
    % always use diagonal 6-9
    bdFlag(4*NT+1:5*NT,2) = bdFlag(t,3);
    bdFlag(5*NT+1:6*NT,4) = bdFlag(t,4);
    bdFlag(6*NT+1:7*NT,3) = bdFlag(t,2);
    bdFlag(7*NT+1:8*NT,1) = bdFlag(t,1);
    % change t in the last
    bdFlag(t,1) = 0;
else
    bdFlag = [];
end
