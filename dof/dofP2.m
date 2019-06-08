function [elem2dof,edge,bdDof] = dofP2(elem)
%% DOFP2 dof structure for P2 element.
%
%  [elem2dof,edge,bdDof] = DOFP2(elem) constructs the dof structure for the
%  quadratic element based on a triangle. elem2dof(t,i) is the global index
%  of the i-th dof of the t-th element. Each triangle contains 6 dofs which
%  are organized according to the order of nodes and edges, i.e.
%  elem2dof(t,1:3) is the pointer to dof on nodes and elem2dof(t,4:6) to
%  three edges. The i-th edge is opposited to the i-th vertex.
%
%  See also dof3P2.
% 
%  Documentation: <a href="matlab:ifem PoissonP2femrate">Quadratic Element
%  for Poisson Equation in 2D</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
[edge, i2, j] = myunique(totalEdge);
N = max(elem(:)); 
NT = size(elem,1);
NE = size(edge,1);
elem2edge = reshape(j,NT,3);
elem2dof = uint32([elem N+elem2edge]);
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
bdEdgeIdx = (i1 == i2);
isBdDof = false(N+2*NE,1);
isBdDof(edge(bdEdgeIdx,:)) = true;   % nodal 
idx = find(bdEdgeIdx);
isBdDof(N+idx) = true;
bdDof = find(isBdDof);