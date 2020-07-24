function M = getmassmatvec(elem2edge,area,Dlambda,elemType,K)
%% GETMASSMATVEC Get the mass matrix of vector finite element space
%
% M = GETMASSMATVEC(elem2edge,area,Dlambda,elemType,K) get mass matrix of
% the finite element space specified by elemType. The coefficient K is
% piecewise constant.
%
% The elemType can be: 
%
% - "RT0": The lowest order Raviart-Thomas element
% - "BDM1": The lowest order Brezzi-Douglas-Michel element
% - "BDM1B": The BDM element enriched by the curl of cubic bubble function
% - "ND0": The lowest order Nedelec element
% - "RT1": Linear-Quadratic Rvairart-Thomas element 
%
% Note that RT0 and ND0 share the same mass matrix.
%
% RT0, BMD1 are created by Ming Wang at July, 2012. Improved by Long Chen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('K','var'), K = []; end
NE = double(max(elem2edge(:)));
NT = size(elem2edge,1);
DiDj = zeros(NT,3,3);
for i = 1:3
    for j = i:3
        if isempty(K)
            DiDj(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2);
        else
            switch size(K,2)
                case 1 % scalar K
                    DiDj(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2)./K;
                case 2 % diagonal matrix
                    DiDj(:,i,j) =  Dlambda(:,1,i).*Dlambda(:,1,j)./K(:,1) ...
                                 + Dlambda(:,2,i).*Dlambda(:,2,j)./K(:,2);
                case 3 % K is 2 by 2 SPD matrix
                    detK = K(:,1).*K(:,2) - K(:,3).^2;
                    DiDj(:,i,j) = (Dlambda(:,1,i).*Dlambda(:,1,j).*K(:,2) ...
                                 + Dlambda(:,2,i).*Dlambda(:,2,j).*K(:,1) ...
                                 - Dlambda(:,1,i).*Dlambda(:,2,j).*K(:,3) ...
                                 - Dlambda(:,2,i).*Dlambda(:,1,j).*K(:,3))./detK;
                    
            end
        end
        DiDj(:,j,i) = DiDj(:,i,j);
    end
end
localEdge = [2 3; 1 3; 1 2]; % ascend ordering

%% RT0 and ND0
if strcmp(elemType,'RT0') || strcmp(elemType,'ND0')
    M = sparse(NE,NE);
    for i = 1:3
        for j = i:3
            % local to global index map and its sign
            ii = double(elem2edge(:,i));
            jj = double(elem2edge(:,j));
            i1 = localEdge(i,1); i2 = localEdge(i,2); % [i1,i2] is the edge opposite to vertex i.
            j1 = localEdge(j,1); j2 = localEdge(j,2);
            % computation of mass matrix --- (phi_i, phi_j)
            Mij = 1/12*area.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                             - (1+(i1==j2))*DiDj(:,i2,j1) ...
                             - (1+(i2==j1))*DiDj(:,i1,j2) ...
                             + (1+(i2==j2))*DiDj(:,i1,j1));
            if (j==i)
                M = M + sparse(ii,jj,Mij,NE,NE);
            else
                M = M + sparse([ii;jj],[jj;ii],[Mij; Mij],NE,NE);
            end
        end
    end
end

%% BDM1
if strcmp(elemType,'BDM1') || strcmp(elemType,'BDM1B') || strcmp(elemType,'RT1')
    M = sparse(2*NE,2*NE);
    for i = 1:3
        for j = i:3
            % local to global index map and its sign
            ii = double(elem2edge(:,i));
            jj = double(elem2edge(:,j));
            i1 = localEdge(i,1); i2 = localEdge(i,2); % [i1,i2] is the edge opposite to vertex i.
            j1 = localEdge(j,1); j2 = localEdge(j,2);
            % computation of mass matrix, note that (rot u, rot v) = (grad u, grad v)
            % (phi_i, phi_j)
            Mij = 1/12*area.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                             - (1+(i1==j2))*DiDj(:,i2,j1) ...
                             - (1+(i2==j1))*DiDj(:,i1,j2) ...
                             + (1+(i2==j2))*DiDj(:,i1,j1));
            if (j==i)
                M = M + sparse(ii,jj,Mij,2*NE,2*NE);
            else
                M = M + sparse([ii;jj],[jj;ii],[Mij; Mij],2*NE,2*NE);
            end
            % (psi_i,psi_j)
            Mij = 1/12*area.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                             + (1+(i1==j2))*DiDj(:,i2,j1) ...
                             + (1+(i2==j1))*DiDj(:,i1,j2) ...
                             + (1+(i2==j2))*DiDj(:,i1,j1));
            if (j==i)
                M = M + sparse(ii+NE,jj+NE,Mij,2*NE,2*NE);
            else
                M = M + sparse([ii;jj]+NE,[jj;ii]+NE,[Mij; Mij],2*NE,2*NE);
            end
        end
    end
    for i = 1:3
        for j = 1:3
            % local to global index map and its sign
            ii = double(elem2edge(:,i));
            jj = double(elem2edge(:,j));
            i1 = localEdge(i,1); i2 = localEdge(i,2);
            j1 = localEdge(j,1); j2 = localEdge(j,2);
            % (psi_i,phi_j)
            Mij = 1/12*area.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                             - (1+(i1==j2))*DiDj(:,i2,j1) ...
                             + (1+(i2==j1))*DiDj(:,i1,j2) ...
                             - (1+(i2==j2))*DiDj(:,i1,j1));
            M = M + sparse([ii+NE;jj],[jj;ii+NE],[Mij; Mij],2*NE,2*NE);
        end
    end
end

%% BDM1B
if strcmp(elemType,'BDM1B')
    newM = sparse(2*NE+NT,2*NE+NT);
    newM(1:2*NE,1:2*NE) = M;
    M = newM;
    % (phi_i, bubble)
    for i = 1:3
        ii = double(elem2edge(:,i)); 
        jj = double((1:NT)');
        i1 = localEdge(i,1); i2 = localEdge(i,2); i3 = 6-i1-i2;
        Mij = 9/10*area.*(...
             dot(Dlambda(:,:,i2),Dlambda(:,:,i3),2)...
            +dot(Dlambda(:,:,i2),Dlambda(:,:,i2),2)...
            -dot(Dlambda(:,:,i1),Dlambda(:,:,i3),2)...
            -dot(Dlambda(:,:,i1),Dlambda(:,:,i1),2));
        if ~isempty(K)
            Mij = Mij./K;
        end        
        M = M + sparse([ii;jj+2*NE],[jj+2*NE;ii],[Mij;Mij],2*NE+NT,2*NE+NT);
    end
    clear ii jj i1 i2 i3 Mij;
    % (psi_i, bubble)
    for i = 1:3
        ii = double(elem2edge(:,i)); 
        jj = double((1:NT)');
        i1 = localEdge(i,1); i2 = localEdge(i,2); i3 = 6-i1-i2;
        Mij = 9/10*area.*(...
             dot(Dlambda(:,:,i2),Dlambda(:,:,i3),2)...
            +dot(Dlambda(:,:,i1),Dlambda(:,:,i3),2)...
            +dot(Dlambda(:,:,i2),Dlambda(:,:,i2),2)...
            +dot(Dlambda(:,:,i1),Dlambda(:,:,i1),2)...
            +dot(Dlambda(:,:,i1),Dlambda(:,:,i2),2));
        if ~isempty(K)
            Mij = Mij./K;
        end        
        M = M + sparse([ii+NE;jj+2*NE],[jj+2*NE;ii+NE],[Mij;Mij],2*NE+NT,2*NE+NT);
    end
    clear ii jj i1 i2 i3 Mij;
    % (bubble, bubble)
    ii = 1+2*NE:NT+2*NE;
    Mij = 81/10*area.*(...
          dot(Dlambda(:,:,1),Dlambda(:,:,1),2)...
         +dot(Dlambda(:,:,1),Dlambda(:,:,2),2)...
         +dot(Dlambda(:,:,1),Dlambda(:,:,3),2)...
         +dot(Dlambda(:,:,2),Dlambda(:,:,2),2)...
         +dot(Dlambda(:,:,2),Dlambda(:,:,3),2)...
         +dot(Dlambda(:,:,3),Dlambda(:,:,3),2));
    if ~isempty(K)
        Mij = Mij./K;
    end        
    M = M + sparse(ii,ii,Mij,2*NE+NT,2*NE+NT);
end

%% RT1
elem2dof = [elem2edge elem2edge+NE (1:NT)'+2*NE (1:NT)'+2*NE+NT];
locBasesIdx = [1 2 0; 1 3 0; 2 3 0; ... % phi and psi
               2 1 3; 3 1 2];           % chi
if strcmp(elemType,'RT1')
    newM = sparse(2*NE+2*NT,2*NE+2*NT);
    newM(1:2*NE,1:2*NE) = M;
    M = newM;
    % (phi_i, chi_j)
    for i = 1:3
        for j = 1:2
            ii = double(elem2dof(:,i));
            jj = double(elem2dof(:,j+6));
            i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); %i3 = locBasesIdx(i,3);
            j1 = locBasesIdx(j,1); j2 = locBasesIdx(j,2); j3 = locBasesIdx(j,3);
            Mij = intlambda([j1,i1,j2],2)*DiDj(:,i2,j3) ...
                 -intlambda([j1,i1,j3],2)*DiDj(:,i2,j2) ...
                 -intlambda([j1,i2,j2],2)*DiDj(:,i1,j3) ...
                 +intlambda([j1,i2,j3],2)*DiDj(:,i1,j2);                        
            M = M + sparse([ii;jj],[jj;ii],[Mij;Mij],2*NE+2*NT,2*NE+2*NT);
        end
    end
    clear ii jj i1 i2 i3 Mij;
    % (psi_i, chi_j)
    for i = 1:3
        for j = 1:2
            ii = double(elem2dof(:,i+3));
            jj = double(elem2dof(:,j+6));
            i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); %i3 = locBasesIdx(i,3);
            j1 = locBasesIdx(j,1); j2 = locBasesIdx(j,2); j3 = locBasesIdx(j,3);
            Mij = intlambda([j1,i1,j2],2)*DiDj(:,i2,j3) ...
                 -intlambda([j1,i1,j3],2)*DiDj(:,i2,j2) ...
                 +intlambda([j1,i2,j2],2)*DiDj(:,i1,j3) ...
                 -intlambda([j1,i2,j3],2)*DiDj(:,i1,j2);                        
            M = M + sparse([ii;jj],[jj;ii],[Mij;Mij],2*NE+2*NT,2*NE+2*NT);
        end
    end
    clear ii jj i1 i2 i3 Mij;
    for i = 1:2
        for j = 1:2
            % (chi_j,chi_i)
            ii = double(elem2dof(:,i+3));
            jj = double(elem2dof(:,j+6));
            i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); i3 = locBasesIdx(i,3);
            j1 = locBasesIdx(j,1); j2 = locBasesIdx(j,2); j3 = locBasesIdx(j,3);
            Mij = intlambda([i1,j1,i2,j2],2)*DiDj(:,i3,j3) ...
                 -intlambda([i1,j1,i2,j3],2)*DiDj(:,i3,j2) ...
                 -intlambda([i1,j1,i3,j2],2)*DiDj(:,i2,j3) ...
                 +intlambda([i1,j1,i3,j3],2)*DiDj(:,i2,j2);                        
            M = M + sparse([ii;jj],[jj;ii],[Mij;Mij],2*NE+2*NT,2*NE+2*NT);
        end
    end
end
