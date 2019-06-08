function [elem2dof,edge] = dof3P2(elem)
%% DOF3P2 dof structure for P2 element in 3-D.
%
% [elem2dof,edge,elem2edge] = DOF3P2(elem) constructs the dof structure for
% the quadratic element based on a tetrahedron. elem2dof(t,i) is the global
% index of the i-th dof of the t-th element. Each tetrahedron contains 10
% dofs which are organized according to the order of nodes and edges, i.e.
% elem2dof(t,1:4) is the pointer to dof on nodes and elem2dof(t,5:10) to
% six edges. The local edges are sorted lexigraphically.
%
% See also dofP2.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

N = max(elem(:)); 
NT = size(elem,1);
totalEdge = uint32([elem(:,[1 2]); elem(:,[1 3]); elem(:,[1 4]); ...
                    elem(:,[2 3]); elem(:,[2 4]); elem(:,[3 4])]);
totalEdge = sort(totalEdge,2);
[edge,i2,j] = myunique(totalEdge);
elem2edge = reshape(j,NT,6);
elem2dof = uint32([elem N+elem2edge]);