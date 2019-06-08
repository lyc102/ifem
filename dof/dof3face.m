function [elem2face,face,elem2faceSign] = dof3face(elem)
%% DOF3FACE dof sturctrue for faces in 3-D.
%
% [elem2face,face,elem2faceSign] = DOF3FACE(elem) constructs data structure
% associated to faces in 3-D. elem is the connectivity matrix for a 3-D
% triangulation. elem2face is the elementwise pointer from elem to face
% indices. face is the face matrix using the ascend ordering. elem2faceSign
% records the consistency of the local faces with induced orientation and
% the global face with ascend orientation.
%
% DOF3FACE is used for elements using dof associated to faces like CR and
% WG elements in 3-D. For elements with orientation such as RT and BDM
% elements, when we use ascend ordering system, no elem2facSign is needed;
% see sc3doc for details. 
%
%  See also dofedge, dof3edge
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

NT = size(elem,1);
totalFace = int32([elem(:,[2 3 4]); elem(:,[1 4 3]); ...
                   elem(:,[1 2 4]); elem(:,[1 3 2])]); % induced ordering
[face, tempvar, j] = myunique(sort(totalFace,2));               
elem2face = uint32(reshape(j,NT,4));
elem2faceSign = int8(reshape(sum(sign(diff(totalFace(:,[1:3,1]),1,2)),2),NT,4));      