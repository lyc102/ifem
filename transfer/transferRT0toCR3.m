function [pro, res] = transferRT0toCR3(node,elem,fixedface)
% Reference: https://arxiv.org/abs/2102.03396
%
%
% See also transferRT0toCR

% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Data Structure and Basis for RT0
[elem2face,face,dofSign] = dof3face(elem);
Dlambda = gradbasis3(node,elem);
NF = size(face,1);
locFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];

%% Assembling
T = sparse(3*NF,NF);
for i = 1:4
    for j = 1:4
        ii = double(elem2face(:,i));
		jj = double(elem2face(:,j));
        j1 = locFace(j,1); j2 = locFace(j,2);  j3 = locFace(j,3);
        % [j1,j2] is the edge opposite to vertex j.
%         sT = (1/2*(i~=j1)*Clambda(:,:,j2) - 1/2*(i~=j2)*Clambda(:,:,j1))...
%              .*repmat(double(dofSign(:,j)),1,2);      % phi_j evaluates at m_i.
        sT = -(Dlambda(:,:,j1) + Dlambda(:,:,j2) + Dlambda(:,:,j3))...
             .*repmat(double(dofSign(:,j)),1,3)/3;      % phi_j evaluates at c_t.
        sT = 1/2*sT;
        T = T + sparse([ii NF+ii 2*NF+ii], [jj jj jj], sT(:), 3*NF, NF);
    end
end
pro = T;
pro([fixedface;fixedface;fixedface],:) = 2*pro([fixedface;fixedface;fixedface],:);
res = pro';
end