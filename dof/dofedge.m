function [elem2edge,edge,elem2edgeSign,edgeSign] = dofedge(elem)
%% DOFEDGE dof structure for edges.
%
% [elem2edge,edge,elem2edgeSign,edgeSign] = DOFEDGE(elem) constructs data
% structure for finite elements associated to edges including CR
% nonconforming element, Ravairt-Thomas element, and Nedelec element etc.
%
% In the input elem is the connectivity matrix for a 2-D triangulation. In
% the output
%
% - elem2edge: the elementwise pointer from elem to edge. In each triangle,
% the opposite indexing is used for its three edges, e.g., elem2edge(t,1)
% is the global index of the first edge consisting of [2 3] vertices of t.
%
% - edge: the edge matrix satisfying edge(:,1)<edge(:,2). The orientation
% of edge is induced by the ordering of vertices, i.e., from the vertex
% with a smaller global index to a bigger one.
%
% - elem2edgeSign: records the consistency of the local and global edge
% orientation. The orientation of local edges is the induced oritentation
% of triangles, i.e., three local edges are: [2 3; 3 1; 1 2]. 
%
% When both elem and local edges are ascend ordered (local edges [2 3; 1 3;
% 1 2] ), elem2edgeSign = [1 1 1]. In this case, elem2edgeSign from DOFEDGE
% is irrelevant and useless. When elem is ascend ordered and local edge is
% induced ordering [2 3; 3 1; 1 2], then elem2edgeSign = [1 -1 1], i.e.,
% only the edge [3 1] is inconsistent.
%
% - edgeSign: As the lcoal edge is the induced ordering, one interior edge
% will be shared by two triangles with opposite orientation. edgeSign equals
% 1 if the edge is consistent with the local edge of its first element,
% equals -1 otherwise.
%
% See also dof3edge
%
% Doc: <a href="matlab:ifem dofedgedoc">dofedgedoc</a>
%
% Modified by Long Chen and Ming Wang.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Generate edge and elem2edge
T = auxstructure(elem);
elem2edge = T.elem2edge;
edge = T.edge;
edge2elem = T.edge2elem;
NT = size(elem,1); 
NE = size(edge,1);

%% The sign of edges
elem2edgeSign = ones(NT,3);
totalEdge = uint32([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
idx = (totalEdge(:,1)>totalEdge(:,2));
elem2edgeSign(idx) = -1;
% isbdEdge = (edge2elem(:,1) == edge2elem(:,2));
% elem2edgeSign(sub2ind([NT,3],edge2elem(:,1),edge2elem(:,3))) = edgeSign;
% idx = sub2ind([NT,3],edge2elem(~isbdEdge,2),edge2elem(~isbdEdge,4));
% elem2edgeSign(idx) = -edgeSign(~isbdEdge);
% Use sub2ind to change the subscript to linear indexing.


%% Consistency of oritentation of edges
% edgeSign equals 1 if the edge is consistent with the local edge of its
% first element, equals -1 otherwise.
if nargout >3
    edgeSign = -ones(NE,1); 
    isConsistentEdge = (elem(sub2ind([NT,3],edge2elem(:,1),mod(edge2elem(:,3),3)+1))== edge(:,1)); 
    edgeSign(isConsistentEdge) = 1;
end
% t = edge2elem(e,1) % t is the first triangle containing e
% k = edge2elem(e,3) % the k-th edge of t is e
% elem(edge2elem(:,1) + NT*mod(edge2elem(:,3),3)) % the first node of
% the k-th of t treating elem as a 1-d array. 
% elem(edge2elem(:,1) + NT*(mod(edge2elem(:,3)+1,3)) % the second node 
% Use sub2ind to change the subscript to linear indexing.
