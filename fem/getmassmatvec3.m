function M = getmassmatvec3(elem2dof,volume,Dlambda,elemType,K)
%% GETMASSMATVEC3 get a mass matrix of finite element spaces in 3D
%
% M = GETMASSMATVEC3(elem2dof,volume,Dlambda,elemType,K) get mass matrix of
% the finite element space specified by elemType. The coefficient K is
% piecewise constant.
%
% The elemType can be: 
%
% - "RT0": The lowest order Raviart-Thomas element
% - "ND0": The lowest order Nedelec element
% - "ND1": The full linear Nedelec element
% - "ND2m": The incomplete quadratic Nedelec element (P\Lambda_2^{-})
% - "ND2": The full quadratic Nedelec element (P\Lambda_2)
%
% Note that for H(div) element, the coefficient is Mij/K while for H(curl)
% element, it is Mij*K. 
%
% See also Maxwell, Maxwell1, Maxwell2, Poisson3RT0
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('K','var'), K = []; end
NT = size(elem2dof,1);

if contains(elemType, 'ND')
    DiDj = zeros(NT,4,4);
    for i = 1:4
        for j = i:4        
            DiDj(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2);
            DiDj(:,j,i) = DiDj(:,i,j);
        end
    end
end

%% ND0: the lowest order edge element
if strcmp(elemType,'ND0')
    NE = double(max(elem2dof(:)));
    locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    ii = zeros(21*NT,1); jj = zeros(21*NT,1); sM = zeros(21*NT,1);
    index = 0;
    for i = 1:6
        for j = i:6
            ii(index+1:index+NT) = double(elem2dof(:,i)); 
            jj(index+1:index+NT) = double(elem2dof(:,j));
            % mass matrix
            % locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            i1 = locEdge(i,1); i2 = locEdge(i,2);
            j1 = locEdge(j,1); j2 = locEdge(j,2);
            Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                               - (1+(i1==j2))*DiDj(:,i2,j1) ...
                               - (1+(i2==j1))*DiDj(:,i1,j2) ...
                               + (1+(i2==j2))*DiDj(:,i1,j1));
            if ~isempty(K)
                Mij = Mij.*K;
            end
            sM(index+1:index+NT) = Mij;
            index = index + NT;
        end
    end
    diagIdx = (ii == jj);   upperIdx = ~diagIdx;
    M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),NE,NE);
    MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),NE,NE);
    M = M + MU + MU';
end

%% ND1: The full linear Nedelec element
if strcmp(elemType,'ND1') 
    NE = double(max(elem2dof(:)));
    Ndof = 2*NE;
    locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    ii = zeros(21*NT,1); jj = zeros(21*NT,1); 
    sMphi = zeros(21*NT,1); sMpsi = zeros(21*NT,1);
    index = 0;
    % two block diagonal matrix for phi and psi
    for i = 1:6
        for j = i:6
            % (phi_i,phi_j)
            ii(index+1:index+NT) = double(elem2dof(:,i)); 
            jj(index+1:index+NT) = double(elem2dof(:,j));
            % mass matrix
            % locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
            i1 = locEdge(i,1); i2 = locEdge(i,2);
            j1 = locEdge(j,1); j2 = locEdge(j,2);
            Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                               - (1+(i1==j2))*DiDj(:,i2,j1) ...
                               - (1+(i2==j1))*DiDj(:,i1,j2) ...
                               + (1+(i2==j2))*DiDj(:,i1,j1));
            if ~isempty(K)
                Mij = Mij.*K;
            end
            sMphi(index+1:index+NT) = Mij;
            % (psi_i,psi_j)
            Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                               + (1+(i1==j2))*DiDj(:,i2,j1) ...
                               + (1+(i2==j1))*DiDj(:,i1,j2) ...
                               + (1+(i2==j2))*DiDj(:,i1,j1));
            if ~isempty(K)
                Mij = Mij.*K;
            end
            sMpsi(index+1:index+NT) = Mij;
            index = index + NT;
        end
    end
    diagIdx = (ii == jj);   upperIdx = ~diagIdx;
    M = sparse([ii(diagIdx); ii(diagIdx)+NE], [jj(diagIdx); jj(diagIdx)+NE],...
               [sMphi(diagIdx); sMpsi(diagIdx)],Ndof,Ndof);
    MU = sparse([ii(upperIdx); ii(upperIdx)+NE], [jj(upperIdx); jj(upperIdx)+NE],...
               [sMphi(upperIdx); sMpsi(upperIdx)],Ndof,Ndof);
    M = M + MU + MU';
    % off-diagonal matrix (psi, phi)
    ii = zeros(36*NT,1); jj = zeros(36*NT,1); ss = zeros(36*NT,1);
    index = 0;
    for i = 1:6
        for j = 1:6
            % local to global index map and its sign
            i1 = locEdge(i,1); i2 = locEdge(i,2);
            j1 = locEdge(j,1); j2 = locEdge(j,2);
            Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                       - (1+(i1==j2))*DiDj(:,i2,j1) ...
                       + (1+(i2==j1))*DiDj(:,i1,j2) ...
                       - (1+(i2==j2))*DiDj(:,i1,j1));
            if ~isempty(K)
                Mij = Mij.*K;
            end
            ii(index+1:index+NT) = double(elem2dof(:,i))+NE; 
            jj(index+1:index+NT) = double(elem2dof(:,j));
            ss(index+1:index+NT) = Mij;
            index = index + NT;
        end
    end
    ML = sparse(ii,jj,ss,Ndof,Ndof);
    M  = M  + ML + ML';    
end

%% ND2m: the incomplete quadratic Nedelec element
% see Maxwell2
% elem2dof = [elem2edge elem2edge+NE elem2face+2*NE elem2face+2*NE+NF];

if strcmp(elemType,'ND2m')
    NE = double(max(elem2dof(:, 1:6), [], "all"));
    NF = double(max(elem2dof(:, 12:16)-2*NE, [], "all"));
    Ndof = 2*(NE+NF);
    locBasesIdx = [1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % Nd0
               1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % grad(P2)
               3 2 4; 3 1 4; 2 1 4; 2 1 3; ...
               4 2 3; 4 1 3; 4 1 2; 3 1 2]; % face bubbles
    ii = zeros(210*NT,1); jj = zeros(210*NT,1);
    sM = zeros(210*NT,1);
    index = 0;
    for i = 1:20
        for j = i:20
            ii(index+1:index+NT) = double(elem2dof(:,i));
            jj(index+1:index+NT) = double(elem2dof(:,j));
            i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); i3 = locBasesIdx(i,3);
            j1 = locBasesIdx(j,1); j2 = locBasesIdx(j,2); j3 = locBasesIdx(j,3);
            Mij = zeros(NT,1);
            if (i<=6) && (j<=6)
                % block 1: (phi_i,phi_j)
                Mij = 1/20*((1+(i1==j1))*DiDj(:,i2,j2) ...
                    - (1+(i1==j2))*DiDj(:,i2,j1) ...
                    - (1+(i2==j1))*DiDj(:,i1,j2) ...
                    + (1+(i2==j2))*DiDj(:,i1,j1));
            end
            if (i<=6) && (7<=j) && (j<=12)
                % block 2: (psi_j,phi_i)
                Mij = 1/20*( (1+(j1==i1))*DiDj(:,j2,i2) ...
                    - (1+(j1==i2))*DiDj(:,j2,i1) ...
                    + (1+(j2==i1))*DiDj(:,j1,i2) ...
                    - (1+(j2==i2))*DiDj(:,j1,i1));

            end
            if (7<=i) && (i<=12) && (7<=j) && (j<=12)
                % block 3: (psi_j,psi_i)
                Mij = 1/20*((1+(i1==j1))*DiDj(:,i2,j2) ...
                    + (1+(i1==j2))*DiDj(:,i2,j1) ...
                    + (1+(i2==j1))*DiDj(:,i1,j2) ...
                    + (1+(i2==j2))*DiDj(:,i1,j1));
            end
            if (i<=6) && (j>12)
                % block 4: (chi_j,phi_i)
                Mij = intlambda([j1,i1,j2],3)*DiDj(:,i2,j3) ...
                    -intlambda([j1,i1,j3],3)*DiDj(:,i2,j2) ...
                    -intlambda([j1,i2,j2],3)*DiDj(:,i1,j3) ...
                    +intlambda([j1,i2,j3],3)*DiDj(:,i1,j2);
            end
            if (7<=i) && (i<=12) && (j>12)
                % block 5: (chi_j,psi_i)
                Mij = intlambda([j1,i1,j2],3)*DiDj(:,i2,j3) ...
                    -intlambda([j1,i1,j3],3)*DiDj(:,i2,j2) ...
                    +intlambda([j1,i2,j2],3)*DiDj(:,i1,j3) ...
                    -intlambda([j1,i2,j3],3)*DiDj(:,i1,j2);
            end
            if (i>12) && (j>12)
                % block 6: (chi_j,chi_i)
                Mij = intlambda([i1,j1,i2,j2],3)*DiDj(:,i3,j3) ...
                    -intlambda([i1,j1,i2,j3],3)*DiDj(:,i3,j2) ...
                    -intlambda([i1,j1,i3,j2],3)*DiDj(:,i2,j3) ...
                    +intlambda([i1,j1,i3,j3],3)*DiDj(:,i2,j2);
            end
            if ~isempty(K)
                Mij = Mij.*K;
            end
            Mij = Mij.*volume;
            sM(index+1:index+NT) = Mij;
            index = index + NT;
        end
    end

    diagIdx = (ii == jj);   upperIdx = ~diagIdx;
    M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
    MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
    M = M + MU + MU';
end

%% ND2: the complete quadratic Nedelec element
% elem2dof = [elem2edge elem2edge+NE elem2edge+2*NE ...
%             elem2face+3*NE elem2face+3*NE+NF elem2face+3*NE+2*NF]
% index 1 to 18: edge basis
% index 19 to 30: face basis
% phi_1 = lambda_i^2 \nabla lambda_j
% phi_2 = lambda_j^2 \nabla lambda_i
% phi_3 = lambda_i lambda_j \nabla (lambda_i - lambda_j)
% psi_{F_{ijk}; 1,2,3} = ijk \in cyc {\lambda_i \lambda_j \nabla \lambda_k}

if strcmp(elemType,'ND2')
    NE = double(max(elem2dof(:, 1:6), [], "all"));
    NF = double(max(elem2dof(:, 19:22)-3*NE, [], "all"));
    Ndof = 3*(NE+NF);
    locBasesIdx = [1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % phi_1
               2 1 0; 3 1 0; 4 1 0; 3 2 0; 4 2 0; 4 3 0; ... % psi_2
               1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % psi_3
               2 3 4; 3 4 2; 4 2 3; ... % face 1
               1 3 4; 3 4 1; 4 1 3; ... % face 2
               1 2 4; 2 4 1; 4 1 2; ... % face 3
               1 2 3; 2 3 1; 3 1 2]; % face 4
    ii = zeros(465*NT,1); jj = zeros(465*NT,1);
    sM = zeros(465*NT,1);
    index = 0;
    for i = 1:30
        for j = i:30
            ii(index+1:index+NT) = double(elem2dof(:,i));
            jj(index+1:index+NT) = double(elem2dof(:,j));
            i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); i3 = locBasesIdx(i,3);
            j1 = locBasesIdx(j,1); j2 = locBasesIdx(j,2); j3 = locBasesIdx(j,3);
            Mij = zeros(NT,1);
            if (i<=12) && (j<=12)
                % block 1: (phi_{1,2},phi_{1,2})
                Mij = intlambda([i1,i1,j1,j1],3)*DiDj(:,i2,j2);
            end
            if (i<=12) && (12<j) && (j<=18)
                % block 2: (phi_{1,2},phi_{3})
                Mij = intlambda([i1,i1,j1,j2],3)*...
                    (DiDj(:,i2,j1) - DiDj(:,i2,j2));
            end
            if (12<i) && (i<=18) && (12<j) && (j<=18)
                % block 3: (phi_{3},phi_{3})
                Mij = intlambda([i1,i2,j1,j2],3)* ...
                    (DiDj(:,i1,j1) + DiDj(:,i2,j2) ...
                    -DiDj(:,i1,j2) - DiDj(:,i2,j1));
            end
            if (i<=12) && (j>18)
                % block 4: (phi_{1,2}, psi_F)
                Mij = intlambda([i1,i1,j1,j2],3)*DiDj(:,i2,j3);
            end
            if (12<i) && (i<=18) && (j>18)
                % block 5: (phi_{3}, psi_F)
                Mij = intlambda([i1,i2,j1,j2],3)*...
                    (DiDj(:,i1,j3) - DiDj(:,i2,j3));
            end
            if (i>18) && (j>18)
                % block 6: (psi_F,psi_F)
                Mij = intlambda([i1,i2,j1,j2],3)*DiDj(:,i3,j3);
            end
            if ~isempty(K)
                Mij = Mij.*K;
            end
            Mij = Mij.*volume;
            sM(index+1:index+NT) = Mij;
            index = index + NT;
        end
    end

    diagIdx = (ii == jj);   upperIdx = ~diagIdx;
    M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
    MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
    M = M + MU + MU';
end

%% RT0: the lowest order face element
if strcmp(elemType,'RT0')
    NF = double(max(elem2dof(:)));
    localFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3]; % ascend ordering
    M = sparse(NF,NF);
    for i = 1:4
        for j = i:4 
            % local to global index map
            ii = double(elem2dof(:,i));
            jj = double(elem2dof(:,j));
            i1 = localFace(i,1); i2 = localFace(i,2); i3 = localFace(i,3);
            j1 = localFace(j,1); j2 = localFace(j,2); j3 = localFace(j,3);
            % computation of mass matrix --- (phi_i, phi_j) 
            Mij = 1/5*volume.*( ...
                  (1+(i1==j1))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
                                   mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
                 +(1+(i1==j2))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
                                   mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
                 +(1+(i1==j3))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
                                   mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)...
                 +(1+(i2==j1))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
                                   mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
                 +(1+(i2==j2))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
                                   mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
                 +(1+(i2==j3))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
                                   mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)...
                 +(1+(i3==j1))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
                                   mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
                 +(1+(i3==j2))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
                                   mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
                 +(1+(i3==j3))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
                                   mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)); 
            if ~isempty(K)
                Mij = Mij./K;
            end
            if (j==i)
                M = M + sparse(ii,jj,Mij,NF,NF);
            else
                M = M + sparse([ii;jj],[jj;ii],[Mij; Mij],NF,NF);        
            end        
        end
    end
end
