function [node,elem,bdEdge,HB] = smeshuniformrefine(node,elem,bdEdge)
%% UNIFORMREFINE refine 2-D triangulation uniformly
%
% See also uniformrefine3, uniformbisect, bisect, bisect3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Construct data structure

T = auxstructure(elem);
edge = T.edge;
elem2edge = T.elem2edge;
edge2elem = T.edge2elem;

N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Add new nodes: for surface mesh, we try to use loop refine to update the
%  location of the new point and old point;

Idxi = elem(edge2elem(:,1)+(edge2elem(:,3)-1)*NT); 
Idxj = elem(edge2elem(:,2)+(edge2elem(:,4)-1)*NT);
node(N+1:N+NE,:) = (3*node(edge(:,1),:) + 3*node(edge(:,2),:) + node(Idxi,:) +node(Idxj,:))/8;
valence = accumarray(elem(:),1,[N,1]);
a = (5/8 - (3/8 + 0.25*cos(2*pi./valence).^2))./valence;
p2p = sparse(double(edge(:,[1,2])),double(edge(:,[2,1])),1, N, N);
node(1:N,:) = repmat((1-valence.*a),1,3).*node(1:N,:) + repmat(a,1,3).*(p2p*node(1:N,:));
HB(:,[1 2 3]) = [(N+1:N+NE)', edge(:,1:2)]; 
edge2newNode = uint32((N+1:N+NE)');

%% Refine each triangle into four triangles as follows
%     3
%    / \
%   5 - 4
%  / \ / \
% 1 - 6 - 2
t = 1:NT;
p(t,1:3) = elem(t,1:3);
p(t,4:6) = edge2newNode(elem2edge(t,1:3));
elem(t,:) = [p(t,1), p(t,6), p(t,5)];
elem(NT+1:2*NT,:) = [p(t,6), p(t,2), p(t,4)];
elem(2*NT+1:3*NT,:) = [p(t,5), p(t,4), p(t,3)];
elem(3*NT+1:4*NT,:) = [p(t,4), p(t,5), p(t,6)];

%% Update boundary edges
if (nargin>=3) && (~isempty(bdEdge))
    bdEdge(NT+1:2*NT,[1 3]) = bdEdge(t,[1 3]); 
    bdEdge(2*NT+1:3*NT,[1 2]) = bdEdge(t,[1 2]); 
    bdEdge(3*NT+1:4*NT,1) = 0;
    bdEdge(t,1) = 0;
else
    bdEdge = [];
end