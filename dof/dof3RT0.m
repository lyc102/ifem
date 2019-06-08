function [elem2face,elem2faceSign,face] = dof3RT0(elem)
%% DOF3RT0 dof sturctrue for the lowest order face element in 3-D.
%
% [elem2face,elem2faceSign,face] = dof3RT0(elem) constructs data structure for
% the BDM face element in 3-D. elem is the connectivity matrix for a 3-D
% triangulation. elem2face is the elementwise pointer from elem to dof
% indices. elem2faceSign records the consistency of the local and global face
% orientation. face is the face matrix.
%
%  See also dofRT0, dof3BDM1, dofBDM1
%
% Added by Ming Wang. Modified from dofstructureface by Long Chen. Add more
% explanation comments by Long Chen. Add regenerate face by Long Chen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

elem = uint32(elem);
%% Auxstructure and global index of faces
NT = size(elem,1);
T = auxstructure3(elem);
elem2face = T.elem2face;

%% Regenerate faces with consistent oritentation
% Note that here we don't use T.face, since the way we assign elem2faceSign i.e.
% the orientation of faces is not sorted by the global indices of vertices.
% locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
face = T.face;
idx = (T.face2elem(:,3) == 1);
face(idx,:) = elem(T.face2elem(idx,1), [2 3 4]);
idx = (T.face2elem(:,3) == 2);
face(idx,:) = elem(T.face2elem(idx,1), [1 4 3]);
idx = (T.face2elem(:,3) == 3);
face(idx,:) = elem(T.face2elem(idx,1), [1 2 4]);
idx = (T.face2elem(:,3) == 4);
face(idx,:) = elem(T.face2elem(idx,1), [1 3 2]);
%% Oritentation of faces
% For each face f of a tetrahedron T, label f with 1 if T is the first
% triangle containing f in the unique command (see auxstructure3).
% According to this rule, the sign of boundary faces are always assigned as
% 1, hence it is no need to consider the sign of boundary faces when we
% handle boundary condtion. The sign will determine the orientation of
% faces.

elem2faceSign = -ones(NT,4,'int8');
elem2faceSign(sub2ind([NT,4],T.face2elem(:,1),T.face2elem(:,3))) = 1;