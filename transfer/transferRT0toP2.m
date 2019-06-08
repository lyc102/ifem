function [pro,res] = transferRT0toP2(node,elem,fixedEdge)
%% TRANSFERRT0TOP2 Transfer RT0 function to P2 function
%
%  Created by Ming Wang, with discussion of Luoping Chen, at July, 2012.
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
T = sparse(2*(NE+N),NE);
for i = 1:3
    ie = double(elem2edge(:,i));
    iv = elem(:,i);
    for j = 1:3
        je = double(elem2edge(:,j));
        j1 = locEdge(j,1); j2 = locEdge(j,2); % [j1,j2] is the edge opposite to vertex j.
        sTc = 1/3*(Clambda(:,:,j2)-Clambda(:,:,j1)).*repmat(double(dofSign(:,j)),1,2); % phi_j evaluates at c_t.
        sTe = 1/2*sTc;
        sTv = repmat(1./elem2NTshareVtx(:,i),1,2).*sTc;
        T = T + sparse([N+ie 2*N+NE+ie], [je je], sTe(:), 2*(N+NE), NE); % edge
        T = T + sparse([iv N+NE+iv], [je je], sTv(:),2*(N+NE), NE);      % vertex
    end
end
T([N+fixedEdge 2*N+NE+fixedEdge],:) = 2*T([N+fixedEdge 2*N+NE+fixedEdge],:);
pro = T;
res = pro';
end
