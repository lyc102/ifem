function [elem2dof,dofSign,edge] = dofBDM1(elem)
%% DOFBDM1 dof structure for BDM1 element
%
% [elem2dof,dofSign,edge] = dofBDM1(elem) constructs data structure for
% the BDM face element in 3-D. elem is the connectivity matrix for a 3-D
% triangulation. elem2dof is the elementwise pointer from elem to dof
% indices. dofSign records the consistency of the local and global edge
% orientation. edge is the edge matrix.
%
% Added by Ming Wang.
%
%  See also dof3BDM1, dofRT0, dof3RT0
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[elem2edge,dofSign,edge] = dofRT0(elem);
NT = size(elem,1);  NE = size(edge,1);
totaledge = uint32([elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])]);
[tempvar,i] = sort(totaledge,2); %#ok<*ASGLU>
[tempvar,j]= sort(i,2);
elem2dof = repmat(elem2edge(:),1,2) + uint32((j-1)*NE);
elem2dof = reshape(elem2dof(:),NT,6);
dofSign = repmat(dofSign,1,2);