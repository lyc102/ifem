function elem2dof = dof3expand(elem2dof, n)
% DOF3EXPAND expands the dof indexing from 1:(num of node/edge/face) to 
% 1:(num of dof per node/edge/face)*(num of node/edge/face), and is used in
% combination with other functions for higher order finite element methods.
% 
% elem2Nd2DoF = dof3expand(elem2edge, 3) changes the dof index from the lowest
% Nedelec space (1 per edge) to quadratic Nedelec element (3 per edge).
%
% elem2faceDoF = dof3expand(elem2face, 12) assigns 12 DoFs per face.
%
% See also transferDG3, dofsign3expand
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%%
d = size(elem2dof, 2);
elem2dofNew = cell(d,1);
assert(iscolumn(elem2dof(:,1)))
for i = 1:d
    elem2dofNew{i} = (elem2dof(:,i)-1)*n + (1:n);
end
elem2dof = horzcat(elem2dofNew{:});
end