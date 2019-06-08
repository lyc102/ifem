function [node,elem] = interfaceUniformrefine(node,elem,interfaceEdge,phi)


%% Construct data structure
totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
[edge, ~, j] = unique(totalEdge,'rows');
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
elem2edge = uint32(reshape(j,NT,3));

isInterfaceNode = false(N,1);
isInterfaceNode(interfaceEdge(:)) = true;

%% Add new nodes: middle points of all edges
tmpNode = (node(edge(:,1),:)+node(edge(:,2),:))/2; 
isNewInterfaceNode = isInterfaceNode(edge(:,1)) & isInterfaceNode(edge(:,2));
a = node(edge(isNewInterfaceNode,1),:) - tmpNode(isNewInterfaceNode,:);
b = node(edge(isNewInterfaceNode,2),:) - tmpNode(isNewInterfaceNode,:);
a = tmpNode(isNewInterfaceNode,:) + [-a(:,2), a(:,1)];
b = tmpNode(isNewInterfaceNode,:) + [-b(:,2), b(:,1)];
tmpNode(isNewInterfaceNode,:) = findintersectbisect(phi,a,b);
node(N+1:N+NE,:) = tmpNode;

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


