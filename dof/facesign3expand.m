function elem2faceSign = facesign3expand(elem2faceSign, n)
% FACESIGN3EXPAND expands the face's signs from (NT, 4) to (NT, 4*n)
% where n is the number of dofs associated with each face, this function is
% used when evaluating boundary integral involving local and global normal
% vectors to faces, and is used in combination with other functions 
% for higher order finite element methods.
% 
% elem2faceSign = dof3expand(elem2faceSign, 12) expands the signs to 12 dof per face.
%
% See also transferDG3, dof3expand
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%%
elem2dofSignNew = cell(n,1);
assert(iscolumn(elem2faceSign(:,1)))
for i = 1:n
    elem2dofSignNew{i} = elem2faceSign;
end
elem2faceSign = horzcat(elem2dofSignNew{:});
end