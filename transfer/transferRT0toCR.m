function [pro, res] = transferRT0toCR(node,elem,fixedEdge)
%% TRANSFERRT0TOCR
%  Created by Ming Wang, at July, 2012, with discussion of Luoping Chen. 

%% Data Structure and Basis for RT0
[elem2edge,edge,dofSign] = dofedge(elem);
Clambda = curlbasis(node,elem);
NE = size(edge,1);
locEdge = [2 3; 3 1; 1 2];

%% Assembling
T = sparse(2*NE,NE);
for i = 1:3
    for j = 1:3
        ii = double(elem2edge(:,i));
		jj = double(elem2edge(:,j));
        j1 = locEdge(j,1); j2 = locEdge(j,2); % [j1,j2] is the edge opposite to vertex j.
%         sT = (1/2*(i~=j1)*Clambda(:,:,j2) - 1/2*(i~=j2)*Clambda(:,:,j1))...
%              .*repmat(double(dofSign(:,j)),1,2);      % phi_j evaluates at m_i.
        sT = (1/3*Clambda(:,:,j2) - 1/3*Clambda(:,:,j1))...
             .*repmat(double(dofSign(:,j)),1,2);      % phi_j evaluates at c_t.
        sT = 1/2*sT;
        T = T + sparse([ii NE+ii], [jj jj], sT(:), 2*NE, NE);
    end
end
pro = T;
pro([fixedEdge NE+fixedEdge],:) = 2*pro([fixedEdge NE+fixedEdge],:);
res = pro';
end