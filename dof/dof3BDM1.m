function [elem2dof,dofSign,face] = dof3BDM1(elem)
%% DOF3BDM1 dof structure for BDM1 element in 3-D.
%
% [elem2dof,dofSign,face] = dof3BDM1(elem) constructs data structure for
% the BDM face element in 3-D. elem is the connectivity matrix for a 3-D
% triangulation. elem2dof is the elementwise pointer from elem to dof
% indices. dofSign records the consistency of the local and global face
% orientation. face is the face matrix.
%
%  See also dof3RT0, dofRT0, dofBDM1
%
% Modified by Ming Wang and Long Chen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[elem2face,dofSign,face] = dof3RT0(elem);
NT = size(elem,1);  NF = size(face,1);
totalFace = uint32([elem(:,[2 3 4]); elem(:,[1 4 3]); ...
                    elem(:,[1 2 4]); elem(:,[1 3 2])]);
[tempvar,i] = sort(totalFace,2);
% Each face in an tetrahedral is assigned with 3 dofs, corresponding to the
% vertices of the face. The ordering is associated with the global
% numbering of the vertices, in the sense that the smallest global
% numbering has local number 1, and the largest global numbering has local
% number 3 in the face.
%   elem2dof = zeros(12*NT,1);
%   elem2dof((i(:,1)-1)*4*NT + (1:4*NT)') = elem2face(:);
%   elem2dof((i(:,2)-1)*4*NT + (1:4*NT)') = elem2face(:)+NF;
%   elem2dof((i(:,3)-1)*4*NT + (1:4*NT)') = elem2face(:)+2*NF;
% Other way instead the above four sentences.
[tempvar,j] = sort(i,2); %#ok<*ASGLU>
elem2dof = repmat(elem2face(:),1,3) + uint32((j-1)*NF);
elem2dof = reshape(elem2dof(:),NT,12);
dofSign = repmat(dofSign,1,3);