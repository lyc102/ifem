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
% - "ND0': The lowest order Nedelec element
%

if ~exist('K','var'), K = []; end
NT = size(elem2dof,1);

%% ND0: the lowest order edge element
if strcmp(elemType,'ND0')
    NE = double(max(elem2dof(:)));
    DiDj = zeros(NT,4,4);
    for i = 1:4
        for j = i:4        
            DiDj(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2);
            DiDj(:,j,i) = DiDj(:,i,j);
        end
    end
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
            Mij = 1/20*volume*4.*( ...
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
                Mij = Mij.*K;
            end
            if (j==i)
                M = M + sparse(ii,jj,Mij,NF,NF);
            else
                M = M + sparse([ii;jj],[jj;ii],[Mij; Mij],NF,NF);        
            end        
        end
    end
end
