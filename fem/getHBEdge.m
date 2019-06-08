function [HBEdge,HBEdgeSign] = getHBEdge(node,elem)
%% GETHBEDGE get HB edge topology for uniformrefine grids
%
%  Created by Ming Wang, at March, 2012.
%
%% Data structure.
Tc = auxstructure(elem);
edgec = Tc.edge;
edge2elemc = Tc.edge2elem;
NEc = size(edge2elemc,1); NTc = size(elem,1);
[nodef,elemf] = uniformrefine(node,elem);
Tf = auxstructure(elemf);
edgef = Tf.edge;
elem2edgef = Tf.elem2edge;
[~,~,~,edgeSignc] = dofedge(elem);
[~,~,~,edgeSignf] = dofedge(elemf);

%% ordering for local edge on fine grid.
%  interior edge
%
%                V2*******************************V4
%                *              /  \              *
%               * *            /    \            *
%              *   *          /      \          *
%             *     1        6        4        *
%            *  t3   *      /    t2    \      *
%           *         *    /            \    *
%          *           *  /              \  *
%         *------7-------@--------8--------*
%        * \             /*               *
%       *   \           /  *             *
%      *     \         /    *           *
%     *       3   t1  5      2   t4    *
%    *         \     /        *       *
%   *           \   /          *     *
%  *             \ /             *  *
% V3*******************************V1
% Explanation:
%   For each edge e, edge2elemc will give the incidence of elements sharing
%   edge e, and also the local incidence of e to each element. This will
%   fix the local index for each refined element.
%
% boundary edge
%
%                V2
%                *
%               * *
%              *   *
%             *     1
%            *  t3   *
%           *         *
%          *           *
%         *------7-------@
%        * \             /*
%       *   \           /  *
%      *     \         /    *
%     *       3   t1  5      2
%    *         \     /  tbd   *
%   *           \   /          *
%  *             \ /             *
% V3*******************************V1

%% Label triangles on fine grid.
t1 = 3*NTc + edge2elemc(:,1);
t2 = 3*NTc + edge2elemc(:,2);
t3 =  mod(edge2elemc(:,3)+1,3)*NTc + edge2elemc(:,1);
t4 =  mod(edge2elemc(:,4)+1,3)*NTc + edge2elemc(:,2);
tbd = mod(edge2elemc(:,3),  3)*NTc + edge2elemc(:,1);

%% Assign edge index
isbdEdgec = (edge2elemc(:,1)==edge2elemc(:,2));
HBEdge = zeros(NEc,8); NTf = size(elemf,1);
% t3
HBEdge(:,1) = elem2edgef(sub2ind([NTf,3],t3(:),edge2elemc(:,3)));            % elem2edgef(t3,edge2elemc(e,3))
% t4
HBEdge(:,2) = elem2edgef(sub2ind([NTf,3],t4(:),edge2elemc(:,4)));            % elem2edgef(t4,edge2elemc(e,4))
% t1
HBEdge(:,3) = elem2edgef(sub2ind([NTf,3],t1(:),edge2elemc(:,3)));            % elem2edgef(t1,edge2elemc(e,3))
HBEdge(:,5) = elem2edgef(sub2ind([NTf,3],t1(:),mod(edge2elemc(:,3),3)+1));   % elem2edgef(t1,edge2elemc(e,3)+1)
HBEdge(:,7) = elem2edgef(sub2ind([NTf,3],t1(:),mod(edge2elemc(:,3)+1,3)+1)); % elem2edgef(t1,edge2elemc(e,3)+2)
% t2
HBEdge(:,4) = elem2edgef(sub2ind([NTf,3],t2(:),edge2elemc(:,4)));            % elem2edgef(t2,edge2elemc(e,4))
HBEdge(:,6) = elem2edgef(sub2ind([NTf,3],t2(:),mod(edge2elemc(:,4),3) +1 )); % elem2edgef(t2,edge2elemc(e,4)+1)
HBEdge(:,8) = elem2edgef(sub2ind([NTf,3],t2(:),mod(edge2elemc(:,4)+1,3)+1)); % elem2edgef(t2,edge2elemc(e,4)+2)
% modify for the boundary edges: tbd
HBEdge(isbdEdgec,2) = elem2edgef(sub2ind([NTf,3],tbd(isbdEdgec),edge2elemc(isbdEdgec,3)));  % elem2edgef(tbd,edge2elemc(e,3))

%% Assign HBEdgeSign
HBEdgeSign = ones(NEc,8);
vec  = node(edgec(:,2),:)-node(edgec(:,1),:);
% parallel edges: HBEdgeSign equals to 1 if (e12, ei) > 0, equals to -1 otherwise.
for i = 1:4 % parallel 
    vefi = nodef(edgef(HBEdge(:,i),2),:)-nodef(edgef(HBEdge(:,i),1),:);
    HBEdgeSign(dot(vec,vefi,2) < 0, i) = -1;
end
% non parallel edges.
% t1 or t2 is always the second element of edge 5, 7, or 6, 8.
HBEdgeSign(:,5) =  edgeSignc.*(-edgeSignf(HBEdge(:,5)));
HBEdgeSign(:,6) = -edgeSignc.*(-edgeSignf(HBEdge(:,6)));
HBEdgeSign(:,7) =  edgeSignc.*(-edgeSignf(HBEdge(:,7)));
HBEdgeSign(:,8) = -edgeSignc.*(-edgeSignf(HBEdge(:,8)));
HBEdgeSign(isbdEdgec,[4 6 8]) = 0; % exclude boundary extra contribution