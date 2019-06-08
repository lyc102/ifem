function [pro,res] = transferRT0toP1b(node,elem)
%% TRANSFERRT0toP1b Transfer RT0 function to P1b function
%
%  Created by Ming Wang at Aug., 2012.
%
%% Data Structure
[elem2edge,edge,dofSign] = dofedge(elem);
NE = size(edge,1); NT = size(elem,1); N = size(node,1);
vtx2elem = sparse(repmat((1:NT)',1,3),elem(:),ones(3*NT,1),NT,N);
NTshareVtx = sum(vtx2elem);
elem2NTshareVtx = NTshareVtx(elem); % #elem in the patch sharing vertex.

%% Basis for RT0
Clambda = curlbasis(node,elem);
locEdge = [2 3; 3 1; 1 2];

%% Assembling
T = sparse(2*(NT+N),NE);
for i = 1:3
    iv = elem(:,i);
    for j = 1:3
        je = double(elem2edge(:,j));
        j1 = locEdge(j,1); j2 = locEdge(j,2); % [j1,j2] is the edge opposite to vertex j.
        sTc = 1/3*(Clambda(:,:,j2)-Clambda(:,:,j1)).*repmat(double(dofSign(:,j)),1,2); % phi_j evaluates at c_t.
        sTv = repmat(1./elem2NTshareVtx(:,i),1,2).*sTc;
        T = T + sparse([iv N+NT+iv], [je je], sTv(:),2*(N+NT), NE);      % vertex
    end
end
pro = T;
res = pro';
end