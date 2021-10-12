function [soln,eqn,info] = DivCurl3PDWGP1(node,elem,bdFlag,pde,option,varargin)
%% DivCurl3PDWGP1 div-curl system equation: 
%       primal-dual weak-Galerkin in 3-D to solve the following system.
%       div(epsilon*u)=f  in \Omega, 
%       curl(u) = g       in \Omega, 
%       Dirichlet boundary condition epsilon*u \cdot n =g_n on \Gamma_0, 
%
% The code is vectorized with no loops in element. Data structure follows iFEM tradition.
%
% P1 element to discretize both primal and dual variables
% lambda_h: 4*NT + 3*NF, local 4 nodal basis in bulk + 3 face basis *4 (no BC)
% u_h: 12*NT, local use Nedelec linear (no BC)
% s_h: s_0 + s_b;  4*NT + 3*NF (simply-connected case), s_b = 0 on Gamma_0
% q_h: q_0 + q_b; 12*NT + 6*NF (for H(curl) vector field) q_b\cross n =0
% q_0 uses Nedelec linear, q_b uses Nedelec linear's tangential projection
%
% Reference: A New Numerical Method for Div-Curl Systems with Low Regularity Assumptions
% S. Cao, C. Wang, J. Wang, https://arxiv.org/abs/2101.03466
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end
if ~isfield(option,'quadorder'), option.quadorder = 2; end % bulk
if ~isfield(option,'fquadorder'), option.fquadorder = 2; end % face
if ~isfield(option,'debug'), option.debug=false; end

[lambdaBulk, wBulk] = quadpts3(option.quadorder);
nQuadBulk = size(lambdaBulk,1);

[lambdaF, wF] = quadpts(option.fquadorder);
nQuadF = size(lambdaF,1);

[lambdaAllFace, wAllFace] = quadpts3face(option.fquadorder);
nQuadAllFace = size(lambdaAllFace,1);
nQuadPerFace = nQuadAllFace/4;

if ~exist('bdFlag','var')
    bdFlag = setboundary3(node,elem,'Dirichlet');
end


%% Mesh related geometric quantities

% local face to vertices indices, positive oriented
face3_to_vert = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
face3_to_edge = [6 5 4; 6 3 2; 5 3 1; 4 2 1];

edge_to_vert = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
face2_to_vert = [2 3; 1 3; 1 2];

[elem2face,face] = dof3face(elem);
elem2face = double(elem2face);

NT = size(elem,1);
NF = size(face,1);

% 12 int dof in each tetra now
elem2Vhdof = reshape(1:12*NT, 12, NT)';

face2Mbdof = reshape(1:3*NF, 3, NF)';

face2Wbdof = reshape(1:6*NF, 6, NF)';

% 4 int dof for M0
elem2M0dof = reshape(1:4*NT, 4, NT)';

[Dphi,volume] = gradbasis3(node,elem);

hK = (6*volume).^(1/3);

[curlPhi,volume] = curlu3(node,elem); % Nd1 extra basis is curl free

faceNormal = cross(node(face(:,2),:) - node(face(:,1),:),...
    node(face(:,3),:) - node(face(:,2),:),2);
faceArea = 0.5*sqrt(sum(faceNormal.^2,2));
faceArea2elem = faceArea(elem2face);

center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
    node(elem(:,3),:) + node(elem(:,4),:))/4;

%% permeability tensor
t = cputime;  % record assembling time

if ~isfield(pde,'Eps'), pde.Eps = []; end

if ~isempty(pde.Eps) && isnumeric(pde.Eps)
    Eps2elem = repmat(pde.Eps, [1,1,NT]);   % Eps is a constant matrix
    Eps2elem = permute(Eps2elem,[3,1,2]); % switch the element idx to the 1st dim
end

if ~isempty(pde.Eps) && ~isnumeric(pde.Eps)
    Eps2elem = zeros(NT,3,3);  %#ok<PREALL>
    % Eps2elem is a NTx3x3 array so that Eps2elem(i,:,:) is the tensor at the
    % center of i-th element
    Epstmp = arrayfun(@(rowidx) pde.Eps(center(rowidx,:)), ...
        1:size(center,1), 'UniformOutput',0);
    Eps2elem = cat(3,Epstmp{:});
    Eps2elem = permute(Eps2elem,[3,1,2]); % switch the element idx to the 1st dim
end

%% compute eps \nabla_w phi
EpsDphi = zeros(NT,3,4); 
for i = 1:4
    EpsDphi(:,:,i) = sum(bsxfun(@times, Eps2elem, Dphi(:,:,i)), 2);
end

%% div(eps lambda_i \nabla lambda_j)
DiDj = zeros(NT,4,4);
for i = 1:4
    for j = i:4        
        DiDj(:,i,j) = dot(EpsDphi(:,:,i),Dphi(:,:,j),2);
        DiDj(:,j,i) = DiDj(:,i,j);
    end
end

%% assembling
% DoF order: # of Dofs
% \lambda_h: 4*NT + 3*NF
% u_h: 12*NT
% s_h: s_0 + s_b;  4*NT + 3*NF, s_b = 0 on Gamma_0
% q_h: q_0 + q_b; 12*NT + 6*NF, q_b\cross n =0 on Gamma
% assembling indices following locally within each group

%% equation 1 for u_h and lambda_h

blockUlambda = sparse(4*NT+3*NF+12*NT,4*NT+3*NF+12*NT);  %#ok<NASGU>
% (u_h, \epsilon\nabla_w phi)
% + \sum h_T^{-1} <\lambda_0 - \lambda_b, \varphi_0 - \varphi_b>_{boundary T}

%% (u_h, \nabla_w phi_{face})
% u_h local basis is now linear Nedelec
% lambda_i grad lambda_j - lambda_j grad lambda_i
% lambda_i grad lambda_j + lambda_j grad lambda_i
% size (4*NT+3*NF, 12*NT)
Gface = sparse(3*NF, 12*NT); %(u_h, \epsilon\nabla_w phi_{face})

for j = 1:12 % ND linear basis q_j
    je = j; if j > 6; je = je-6; end
    
    for f = 1:4 % face 1 to 4
        % Cij = (q_j, eps \nabla_w \phi_{F_i})
        % = <\phi_{F_i}, u_h\cdot n_F>
        % n_{F_i} = -3|K| Dphi_{F_i}/|F_i|
        % q_j = lambda_{j1} \nabla \lambda_{j2} - lambda_{j2} \nabla \lambda_{j1}
        % q_j = lambda_{j1} \nabla \lambda_{j2} + lambda_{j2} \nabla \lambda_{j1}
        
        j1 = edge_to_vert(je,1); j2 = edge_to_vert(je,2);
        
        for i = 1:3 % each face has 3 DoFs for scalar
            for p = 1:nQuadF
                pf = nQuadF*(f-1) + p;
                
                if j <= 6
                    uhp = lambdaAllFace(pf,j1)*Dphi(:,:,j2) ...
                        -lambdaAllFace(pf,j2)*Dphi(:,:,j1);
                else
                    uhp = lambdaAllFace(pf,j1)*Dphi(:,:,j2) ...
                        +lambdaAllFace(pf,j2)*Dphi(:,:,j1);
                end
                
                uhpdotn = dot(uhp, -3*EpsDphi(:,:,f), 2).*volume;
                
                Gij = wF(p)*uhpdotn*lambdaF(p,i);
                
                ii = face2Mbdof(elem2face(:,f), i);
                % old: ii = 3*(elem2face(:,f)-1) + i;
                jj = elem2Vhdof(:,j);
                Gface = Gface + sparse(ii,jj,Gij,3*NF,12*NT);
            end
        end
    end
end

%% (u_h, \nabla_w phi_{bulk})
Gbulk = sparse(4*NT, 12*NT); 
%(u_h, \epsilon\nabla_w phi_{bulk}) = -(div (\epsilon u_h),  phi_{bulk})

for j = 1:6 % ND linear basis q_j, the non-div free part
    j1 = edge_to_vert(j,1); j2 = edge_to_vert(j,2);
    for i = 1:4 % each K has 4 DoFs
        divEpsUhPhi0 = -2*DiDj(:, j1, j2).*volume/4;
        
        ii = elem2M0dof(:, i);
        jj = elem2Vhdof(:,j+6);
        Gbulk = Gbulk + sparse(ii,jj,divEpsUhPhi0,4*NT,12*NT);
    end
end

%% (u_h, \nabla_w phi)
Ugradphi = [Gbulk; Gface];

%% stabilization 
% \sum h_T^{-1} <\lambda_0 - \lambda_b, \varphi_0 - \varphi_b>_{boundary T}
Stab_Mh = sparse(4*NT+3*NF,4*NT+3*NF); %#ok<NASGU>

% \sum h_T <s_0 - s_b, r_0 - r_b>_{boundary T}
Stab_Sh = sparse(4*NT+3*NF,4*NT+3*NF);  %#ok<NASGU>

%% <\lambda_0, \lambda_0>
M00 = sparse(4*NT, 4*NT);
S00 = sparse(4*NT, 4*NT);
for i = 1:4
    for j = 1:4
        M00ijface = 0;
        S00ijface = 0;
        for p = 1:nQuadAllFace
            Sij = lambdaAllFace(p, i)*lambdaAllFace(p, j);
            f = ceil(p/nQuadPerFace); % current face
            M00ijface = M00ijface + wAllFace(p)*Sij*faceArea2elem(:,f)./hK;
            S00ijface = S00ijface + wAllFace(p)*Sij*faceArea2elem(:,f).*hK;
        end
        ii = elem2M0dof(:,i);
        jj = elem2M0dof(:,j);
        M00 = M00 + sparse(ii,jj,M00ijface,4*NT,4*NT);
        S00 = S00 + sparse(ii,jj,S00ijface,4*NT,4*NT);
    end
end

%% <\lambda_b, \lambda_0>
Mb0 = sparse(4*NT, 3*NF); %% 
Sb0 = sparse(4*NT, 3*NF);

for f = 1:4 % face 1 to 4
    for j = 1:4 % P1 basis in bulk
        for i = 1:3 % 3 basis functions on each face
            Mb0ijface = 0;
            Sb0ijface = 0;
            for p = 1:nQuadF
                pf = nQuadF*(f-1) + p; 
                Sb0ij = lambdaAllFace(pf,j)*lambdaF(p,i);
                
                Mb0ijface = Mb0ijface + wF(p)*Sb0ij.*faceArea2elem(:,f)./hK;
                Sb0ijface = Sb0ijface + wF(p)*Sb0ij.*faceArea2elem(:,f).*hK;
            end
            ii = face2Mbdof(elem2face(:,f), i);
            % old: ii = 3*(elem2face(:,f)-1) + i;
            % k-th (global face indexing) face has 3 DoFs
            jj = elem2M0dof(:,j);
            Mb0 = Mb0 + sparse(jj,ii,Mb0ijface,4*NT,3*NF);
            Sb0 = Sb0 + sparse(jj,ii,Sb0ijface,4*NT,3*NF);
        end
    end
end

%% <\lambda_b, \lambda_b>
Mbb = sparse(3*NF, 3*NF);
Sbb = sparse(3*NF, 3*NF); 
for f = 1:4 % face 1 to 4
    for i = 1:3
        for j = 1:3
            ii = face2Mbdof(elem2face(:,f), i);
            jj = face2Mbdof(elem2face(:,f), j);
            Lbij = ((i==j)+1)/12*faceArea2elem(:,f);
            Mbb = Mbb + sparse(ii, jj, Lbij./hK, 3*NF, 3*NF);
            Sbb = Sbb + sparse(ii, jj, Lbij.*hK, 3*NF, 3*NF);
        end
    end
end

Stab_Mh = [M00 -Mb0; -Mb0' Mbb];
Stab_Sh = [S00 -Sb0; -Sb0' Sbb];

if option.debug
    % Stablization on M_h should have null space dimension 
    % v_0 = v_b, i.e., linear continuous P1, number of node (check)
    disp(svds(Stab_Sh,size(node,1),'smallest'))
end

%% Equation 1 for u_h and q_h = {q_0, q_b}

%%  (u_h, \nabla_w \times psi_{face})
% u_h basis: linear Nedelec
% psi_{face} basis are linear Nedelec's tangential component
% size: (6*NF, 12*NT)

D = sparse(6*NF,12*NT);
for f = 1:4 % face 1 to 4
    f2v = face3_to_vert(f,:); % local vertices on face F

    for j = 1:12 % ND linear basis q_j
        je = j; if j > 6; je = je-6; end
        j1 = edge_to_vert(je,1); j2 = edge_to_vert(je,2);

        for e = 1:6 % 6 basis vectors on each face
            % Dije = (q_j, \nabla_w \times \psi_{F_i,e})
            % = (q_j, -\psi_{F_i,e}\times n)_{\partial K}
            ee = e; if e > 3; ee = e-3; end
            e1 = face2_to_vert(ee,1); e2 = face2_to_vert(ee,2);
            
            Upsibijface = 0;
            for p = 1:nQuadF
                pf = nQuadF*(f-1) + p;
                
                if j <=6
                    uhp = lambdaAllFace(pf,j1)*Dphi(:,:,j2)-lambdaAllFace(pf,j2)*Dphi(:,:,j1);
                else
                    uhp = lambdaAllFace(pf,j1)*Dphi(:,:,j2)+lambdaAllFace(pf,j2)*Dphi(:,:,j1);
                end
                
                if e <= 3
                    psibf = lambdaF(p,e1)*Dphi(:,:,f2v(e2))-lambdaF(p,e2)*Dphi(:,:,f2v(e1));
                else
                    psibf = lambdaF(p,e1)*Dphi(:,:,f2v(e2))+lambdaF(p,e2)*Dphi(:,:,f2v(e1));
                end

                psibcrossn = mycross(psibf, Dphi(:,:,f));
                
                Upsibijface = Upsibijface + wF(p)*3*dot(uhp,psibcrossn,2).*volume./hK;
                
            end
            ii = face2Wbdof(elem2face(:,f), e);
            % ii = 6*(elem2face(:,f)-1) + e;
            % k-th (global face indexing) face has 3 DoFs
            jj = elem2Vhdof(:,j);
            D = D + sparse(ii,jj,Upsibijface,6*NF,12*NT);
        end
    end
end

%% (u_h, curl_w psi_{bulk}) = (curl u_h, psi_{bulk})
C = sparse(12*NT,12*NT); 

for i = 1:12 % Nd linear basis for psi_bulk
    ie = i; if i > 6; ie = ie-6; end
    i1 = edge_to_vert(ie,1); i2 = edge_to_vert(ie,2);
    
    % lambda_{i1} \nabla \lambda_{i2} - lambda_{i2} \nabla \lambda_{i1}
    for j = 1:6 % Nd0 basis for u_h
        if i <= 6
            Cij = dot(Dphi(:,:,i2)-Dphi(:,:,i1), curlPhi(:,:,j), 2).*volume/4;
        else
            Cij = dot(Dphi(:,:,i2)+Dphi(:,:,i1), curlPhi(:,:,j), 2).*volume/4;
        end
        
        ii = elem2Vhdof(:,i);
        jj = elem2Vhdof(:,j);
        C = C + sparse(ii,jj,Cij,12*NT,12*NT);
    end
end

%% (u_h, curl_w psi): (12*NT+6*NF, 12*NT)
Ucurlpsi = [C; D];

%%  (\psi, \epsilon\nabla_w s_h)
Psigrads = sparse(12*NT+6*NF, 4*NT+3*NF);

%% E = (psi_0, eps \nabla_w s_{face, bulk})
E = sparse(12*NT, 4*NT+3*NF); %#ok<NASGU>

%% (psi_0, eps \nabla_w s_{bulk}) = -(div (\epsilon psi_0),  s_{bulk})
Ebulk = sparse(12*NT, 4*NT);
for j = 1:6 % ND linear basis q_j, the non-div free part
    j1 = edge_to_vert(j,1); j2 = edge_to_vert(j,2);
    for i = 1:4 % each K has 4 DoFs
        divEpsPsi0s0 = -2*DiDj(:, j1, j2).*volume/4;
        
        jj = elem2Vhdof(:,j+6);
        ii = elem2M0dof(:, i);
        Ebulk = Ebulk + sparse(jj,ii,divEpsPsi0s0,12*NT,4*NT);
    end
end


%% (psi_0, eps \nabla_w s_{F_j}) = <\epsilon psi_0 \cdot n,  s_{F_j})>
Eface = sparse(12*NT, 3*NF);
for f = 1:4 % face 1 to 4
    % n_{F_i} = -3|K| Dphi_{F_i}/|F_i|
    normal = -3*EpsDphi(:,:,f).*repmat(volume, 1, 3);
    for j = 1:12 % ND linear basis q_j
        je = j; if j > 6; je = je-6; end
        j1 = edge_to_vert(je,1); j2 = edge_to_vert(je,2);
        
        
        % Cij = (psi_0, eps \nabla_w s_{F_j})
        % = <s_{F_j}, eps u_h\cdot n_{F_j}>
        
        % psi_j = lambda_{j1} \nabla \lambda_{j2} - lambda_{j2} \nabla \lambda_{j1}
        % psi_j = lambda_{j1} \nabla \lambda_{j2} + lambda_{j2} \nabla \lambda_{j1}
        
        for i = 1:3 % each face has 3 DoFs
            spsi0dotnij = 0;
            for p = 1:nQuadF
                pf = nQuadF*(f-1) + p;
                
                if j <= 6
                    psi0p = lambdaAllFace(pf,j1)*Dphi(:,:,j2) ...
                        -lambdaAllFace(pf,j2)*Dphi(:,:,j1);
                else
                    psi0p = lambdaAllFace(pf,j1)*Dphi(:,:,j2) ...
                        +lambdaAllFace(pf,j2)*Dphi(:,:,j1);
                end
                
                psi0pdotn = dot(psi0p, normal, 2);
                
                spsi0dotnij = spsi0dotnij + wF(p)*psi0pdotn*lambdaF(p,i);   
            end
            ii = face2Mbdof(elem2face(:,f), i);
            jj = elem2Vhdof(:,j);
            Eface = Eface + sparse(jj,ii,spsi0dotnij,12*NT,3*NF);
        end
    end
end

%% Eface is Gface's transpose, the code above is for booking-keeping
Eface = Gface';
E = [Ebulk, Eface];
Psigrads(1:12*NT,:) = E;

%% \sum h_T^{-1} <(q_0-q_b)\times n, (\psi_0-\psi_b)\times n>_{boundary T}

%% \sum h_T^{-1} <q_0\times n, \psi_0 \times n>_{boundary T}
% new implementation following face by face integral

WS00 = sparse(12*NT, 12*NT);
% edge_to_vert = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
% face2_to_vert = [2 3; 1 3; 1 2];

for f = 1:4
    nface = -3*repmat(volume,1,3).*Dphi(:,:,f); % scaled exterior n|F|
    loc_f2e = face3_to_edge(f,:);  % locally on each face the edge's indexing
    
    for pp = 1:nQuadF
        for i = 1:3 % q_0
            for j = 1:3 % psi_0
                ii = loc_f2e(i); 
                % locally edge opposite to the i-th vertex on the current face
                jj = loc_f2e(j); 
                i1 = edge_to_vert(ii,1);
                i2 = edge_to_vert(ii,2);
                j1 = edge_to_vert(jj,1);
                j2 = edge_to_vert(jj,2);
                
                vi1 = face2_to_vert(i, 1);
                vi2 = face2_to_vert(i, 2);
                
                vj1 = face2_to_vert(j, 1);
                vj2 = face2_to_vert(j, 2);
                
                phi_i = lambdaF(pp,vi1)*Dphi(:,:,i2) ...
                    - lambdaF(pp,vi2)*Dphi(:,:,i1);
                
                phi_j = lambdaF(pp,vj1)*Dphi(:,:,j2) ...
                    - lambdaF(pp,vj2)*Dphi(:,:,j1);
                
                psi_i = lambdaF(pp,vi1)*Dphi(:,:,i2) ...
                    + lambdaF(pp,vi2)*Dphi(:,:,i1);
                
                psi_j = lambdaF(pp,vj1)*Dphi(:,:,j2) ...
                    + lambdaF(pp,vj2)*Dphi(:,:,j1);
                
                phiicrossn = mycross(phi_i, nface);
                phijcrossn = mycross(phi_j, nface);
                psiicrossn = mycross(psi_i, nface);
                psijcrossn = mycross(psi_j, nface);
                
                phiiphij = wF(pp)*dot(phiicrossn,phijcrossn,2)./faceArea2elem(:,f);
                WS00 = WS00 + sparse(elem2Vhdof(:,ii),elem2Vhdof(:,jj),...
                    phiiphij,12*NT,12*NT);
                
                phiipsij = wF(pp)*dot(phiicrossn,psijcrossn,2)./faceArea2elem(:,f);
                WS00 = WS00 + sparse(elem2Vhdof(:,ii),elem2Vhdof(:,jj+6),...
                    phiipsij,12*NT,12*NT);
                
                psiiphij = wF(pp)*dot(psiicrossn,phijcrossn,2)./faceArea2elem(:,f);
                WS00 = WS00 + sparse(elem2Vhdof(:,ii+6),elem2Vhdof(:,jj),...
                    psiiphij,12*NT,12*NT);
                
                psiipsij = wF(pp)*dot(psiicrossn,psijcrossn,2)./faceArea2elem(:,f);
                WS00 = WS00 + sparse(elem2Vhdof(:,ii+6),elem2Vhdof(:,jj+6),...
                    psiipsij,12*NT,12*NT);
            end
        end
    end
end

%%
% \sum h_T^{-1} <q_b\times n, \psi_0 \times n>_{boundary T}
WSb0 = sparse(12*NT, 6*NF);

% edge_to_vert = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
% face2_to_vert = [2 3; 1 3; 1 2];
% face3_to_edge = [6 5 4; 6 3 2; 5 3 1; 4 2 1];

for f = 1:4
    nface = -3*repmat(volume,1,3).*Dphi(:,:,f); % scaled exterior n|F|
    loc_f2e = face3_to_edge(f,:);  % locally on each face the edge's indexing
    
    for pp = 1:nQuadF
        for i = 1:3 % psi_0
            for j = 1:3 % q_b
                ii = loc_f2e(i); 
                % locally edge opposite to the i-th vertex on the current face
                jj = loc_f2e(j);
                
                i1 = edge_to_vert(ii,1);
                i2 = edge_to_vert(ii,2);
                j1 = edge_to_vert(jj,1);
                j2 = edge_to_vert(jj,2);
                
                vi1 = face2_to_vert(i, 1);
                vi2 = face2_to_vert(i, 2);
                
                vj1 = face2_to_vert(j, 1);
                vj2 = face2_to_vert(j, 2);
                
                psi01i = lambdaF(pp,vi1)*Dphi(:,:,i2) ...
                    - lambdaF(pp,vi2)*Dphi(:,:,i1);
                
                qb1j = lambdaF(pp,vj1)*Dphi(:,:,j2) ...
                    - lambdaF(pp,vj2)*Dphi(:,:,j1);
                
                psi02i = lambdaF(pp,vi1)*Dphi(:,:,i2) ...
                    + lambdaF(pp,vi2)*Dphi(:,:,i1);
                
                qb2j = lambdaF(pp,vj1)*Dphi(:,:,j2) ...
                    + lambdaF(pp,vj2)*Dphi(:,:,j1);
                
                psi01icrossn = mycross(psi01i, nface);
                qb1jcrossn = mycross(qb1j, nface);
                psi02icrossn = mycross(psi02i, nface);
                qb2jcrossn = mycross(qb2j, nface);
                
                qb1j_psi01i = wF(pp)*dot(qb1jcrossn,psi01icrossn,2)./faceArea2elem(:,f);
                WSb0 = WSb0 + sparse(elem2Vhdof(:,ii),...
                                     face2Wbdof(elem2face(:,f), j),...
                                     qb1j_psi01i, 12*NT,6*NF);
                
                qb1j_psi02i = wF(pp)*dot(qb1jcrossn,psi02icrossn,2)./faceArea2elem(:,f);
                WSb0 = WSb0 + sparse(elem2Vhdof(:,ii+6), ...
                                     face2Wbdof(elem2face(:,f), j),...
                                     qb1j_psi02i, 12*NT,6*NF);
                
                qb2j_psi01i = wF(pp)*dot(qb2jcrossn,psi01icrossn,2)./faceArea2elem(:,f);
                WSb0 = WSb0 + sparse(elem2Vhdof(:,ii), ...
                                     face2Wbdof(elem2face(:,f), j+3),...
                                     qb2j_psi01i, 12*NT,6*NF);
                
                qb2j_psi02i = wF(pp)*dot(qb2jcrossn,psi02icrossn,2)./faceArea2elem(:,f);
                WSb0 = WSb0 + sparse(elem2Vhdof(:,ii+6), ...
                                     face2Wbdof(elem2face(:,f), j+3),...
                                     qb2j_psi02i,12*NT,6*NF);
            end
        end
    end
end

%% \sum h_T^{-1} <q_b\times n, \psi_b \times n>_{boundary T}
WSbb = sparse(6*NF, 6*NF);
for f = 1:4
    nface = -3*repmat(volume,1,3).*Dphi(:,:,f); % scaled exterior n|F|
    loc_f2e = face3_to_edge(f,:);  % locally on each face the edge's indexing
    
    for pp = 1:nQuadF
        for i = 1:3 % psi_b
            for j = 1:3 % q_b
                ii = loc_f2e(i); 
                % locally edge opposite to the i-th vertex on the current face
                jj = loc_f2e(j);
                
                i1 = edge_to_vert(ii,1);
                i2 = edge_to_vert(ii,2);
                j1 = edge_to_vert(jj,1);
                j2 = edge_to_vert(jj,2);
                
                vi1 = face2_to_vert(i, 1);
                vi2 = face2_to_vert(i, 2);
                
                vj1 = face2_to_vert(j, 1);
                vj2 = face2_to_vert(j, 2);
                
                psib1i = lambdaF(pp,vi1)*Dphi(:,:,i2) ...
                    - lambdaF(pp,vi2)*Dphi(:,:,i1);
                
                qb1j = lambdaF(pp,vj1)*Dphi(:,:,j2) ...
                    - lambdaF(pp,vj2)*Dphi(:,:,j1);
                
                psib2i = lambdaF(pp,vi1)*Dphi(:,:,i2) ...
                    + lambdaF(pp,vi2)*Dphi(:,:,i1);
                
                qb2j = lambdaF(pp,vj1)*Dphi(:,:,j2) ...
                    + lambdaF(pp,vj2)*Dphi(:,:,j1);
                
                psib1icrossn = mycross(psib1i, nface);
                qb1jcrossn = mycross(qb1j, nface);
                psib2icrossn = mycross(psib2i, nface);
                qb2jcrossn = mycross(qb2j, nface);
                
                qb1j_psib1i = wF(pp)*dot(qb1jcrossn,psib1icrossn,2)./faceArea2elem(:,f);
                WSbb = WSbb + sparse(face2Wbdof(elem2face(:,f), i),...
                                     face2Wbdof(elem2face(:,f), j),...
                                     qb1j_psib1i, 6*NF, 6*NF);
                
                qb1j_psib2i = wF(pp)*dot(qb1jcrossn,psib2icrossn,2)./faceArea2elem(:,f);
                WSbb = WSbb + sparse(face2Wbdof(elem2face(:,f), i+3), ...
                                     face2Wbdof(elem2face(:,f), j),...
                                     qb1j_psib2i, 6*NF, 6*NF);
                
                qb2j_psib1i = wF(pp)*dot(qb2jcrossn,psib1icrossn,2)./faceArea2elem(:,f);
                WSbb = WSbb + sparse(face2Wbdof(elem2face(:,f), i), ...
                                     face2Wbdof(elem2face(:,f), j+3),...
                                     qb2j_psib1i, 6*NF, 6*NF);
                
                qb2j_psib2i = wF(pp)*dot(qb2jcrossn,psib2icrossn,2)./faceArea2elem(:,f);
                WSbb = WSbb + sparse(face2Wbdof(elem2face(:,f), i+3), ...
                                     face2Wbdof(elem2face(:,f), j+3),...
                                     qb2j_psib2i, 6*NF, 6*NF);
            end
        end
    end
end


%%
Stab_Wh = [WS00 -WSb0; -WSb0' WSbb]; 
% Stab_Wh should have null space dimension = continuous Nedelec linear
if option.debug
    [~,edge] = dof3edge(elem);
    NE = size(edge,1);
    fprintf('WS00 symmetric: %d \n', issymmetric(WS00));
    fprintf('WSbb symmetric: %d \n', issymmetric(WSbb));
    fprintf('Linear Nedelec dim: %d \n', NE);
    fprintf('Null space dim of W_h stab: %d \n', size(Stab_Wh,1) - sprank(Stab_Wh));
end

%% overall stiffness matrix for the system
blockUlambda = [Stab_Mh Ugradphi; Ugradphi' sparse(12*NT,12*NT)];
blockOffDiag = sparse(4*NT+3*NF+12*NT+6*NF, 4*NT+3*NF+12*NT);
blockOffDiag(4*NT+3*NF+1:end, 4*NT+3*NF+1:end) = Ucurlpsi;
blockQs = [-Stab_Sh Psigrads'; Psigrads Stab_Wh];
bigA = [blockUlambda blockOffDiag'; blockOffDiag blockQs];

%% verifying kernel of bigA = span{lambda constant, s constant}
% (lambda 4*NT+3*NF, u_h 12*NT, s 4*NT+3*NF, q 12*NT+6*NF)

if option.debug
    kerA = [ones(4*NT+3*NF,1);
        zeros(12*NT,1);
        ones(4*NT+3*NF,1);
        zeros(12*NT+6*NF,1)];
    fprintf('Norm of A*(constant in lambda and s) is %6g \n', norm(bigA*kerA));
    fprintf('Size of A: (%d, %d) \n', size(bigA));
end
%% right hand sides

%% -(f,\varphi_0) P1, same with linear Lagrange

ft = zeros(NT,4);
for p = 1:nQuadBulk
    pxyz = lambdaBulk(p,1)*node(elem(:,1),:) ...
        + lambdaBulk(p,2)*node(elem(:,2),:) ...
        + lambdaBulk(p,3)*node(elem(:,3),:) ...
        + lambdaBulk(p,4)*node(elem(:,4),:);
    fp = pde.f(pxyz);
    for j = 1:4
        ft(:,j) = ft(:,j) + wBulk(p)*lambdaBulk(p,j)*fp;
    end
end
ft = ft.*repmat(volume,1,4);
Fphi0 = -accumarray(elem2M0dof(:),ft(:),[4*NT 1]);

%% \sum_{T} <eps u\cdot n, \varphi_b>_{\partial T} near \partial \Omega
GNphib = zeros(3*NF,1);

for i = 1:4 % 4 faces
    isBdElem = find(bdFlag(:,i) == 1); 
    % if an element has a i-th face as a boundary face 

    if ~isempty(isBdElem)
        bdNormal = -3*repmat(volume(isBdElem,:),1,3).*Dphi(isBdElem,:,:); 
        % scaled normal with face area built-in n|F|
        loc_face = face3_to_vert(i,:);
        
        for p = 1:nQuadF
            pxyz = lambdaF(p,1)*node(elem(isBdElem,loc_face(1)),:) ...
                + lambdaF(p,2)*node(elem(isBdElem,loc_face(2)),:)...
                + lambdaF(p,3)*node(elem(isBdElem,loc_face(3)),:);
            
            u_bd = pde.exactu(pxyz);
            EpsU = sum(bsxfun(@times, Eps2elem(isBdElem,:,:), u_bd), 2);
            EpsUdotN = dot(squeeze(EpsU), bdNormal(:,:,i), 2);
            for j = 1:3
                % jj = 3*(elem2face(isBdElem,i)-1) + j;
                jj = face2Mbdof(elem2face(isBdElem,i), j);
                GNphib = GNphib + accumarray(jj, wF(p)*lambdaF(p,j)*EpsUdotN, [3*NF 1]);
            end 
        end
    end
end

%% (g, \psi_0)
Gpsi0 = zeros(12*NT,1); %#ok<PREALL>
gt = zeros(NT,12);
for p = 1:nQuadBulk
    pxyz = lambdaBulk(p,1)*node(elem(:,1),:) ...
        + lambdaBulk(p,2)*node(elem(:,2),:) ...
        + lambdaBulk(p,3)*node(elem(:,3),:) ...
        + lambdaBulk(p,4)*node(elem(:,4),:);
    ghp = pde.g(pxyz);
    %   loc_e2v = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    for j = 1:12
        je = j; if j > 6; je = je-6; end
        j1 = edge_to_vert(je,1); j2 = edge_to_vert(je,2);
        if j <= 6
            % phi_k = lambda_i Dphi_j - lambda_j Dphi_i;
            phi_j = lambdaBulk(p,j1)*Dphi(:,:,j2)-lambdaBulk(p,j2)*Dphi(:,:,j1);
        else
            phi_j = lambdaBulk(p,j1)*Dphi(:,:,j2)+lambdaBulk(p,j2)*Dphi(:,:,j1);
        end
        rhs = dot(phi_j,ghp,2);
        gt(:,j) = gt(:,j) + wBulk(p)*rhs;
    end
end
gt = gt.*repmat(volume,1,12);
Gpsi0 = accumarray(elem2Vhdof(:),gt(:),[12*NT 1]);

%% RHS
Rhs = [Fphi0; GNphib]; % test: phi_0, phi_b
Rhs = [Rhs; zeros(12*NT, 1)]; % test: v
Rhs = [Rhs; zeros(4*NT, 1)]; % test: r_0
Rhs = [Rhs; zeros(3*NF, 1)]; % test: r_b
Rhs = [Rhs; Gpsi0]; % test: psi_0
Rhs = [Rhs; zeros(6*NF, 1)]; % test: psi_b

%% set up free dof
% DoF order: # of Dofs
% prefix with "idx": global indices
% the DoF variable itself is boolean (fastest on MATLAB)
% lambda_h: 4*NT + 3*NF, local 4+3*4 (no BC)
% u_h: 12*NT (no BC)
% s_h: s_0 + s_b;  4*NT + 3*NF (simply-connected case), s_b = 0 on Gamma_0
% q_h: q_0 + q_b; 12*NT + 6*NF (for H(curl) vector field) q_b\cross n =0

isBdFace = false(NF,1);
isBdFace(elem2face(bdFlag(:) == 1)) = true;
idxBdFace = find(isBdFace);

idxBdFaceLinear = [3*idxBdFace-2 3*idxBdFace-1 3*idxBdFace]';
idxBdFaceLinear = idxBdFaceLinear(:);

idxBdFaceVec = [6*idxBdFace-5 6*idxBdFace-4 6*idxBdFace-3 ...
                6*idxBdFace-2 6*idxBdFace-1 6*idxBdFace]';
idxBdFaceVec = idxBdFaceVec(:);

% ({lambda_0, lambda_b}, u_h, {s_0, s_b}, {q_0, q_b})
NDof = (4*NT + 3*NF) + 12*NT + (4*NT+3*NF) + (12*NT + 6*NF);
freeDof = true(NDof,1);

% boundary DoF indices for s_b
idxBdSb = idxBdFaceLinear + ((4*NT + 3*NF) + 12*NT + 4*NT); 

% boundary DoF indices for q_b
idxBdQb = idxBdFaceVec + ((4*NT + 3*NF) + 12*NT + (4*NT+3*NF) + 12*NT); 

% since for lambda_0 and s_0, adding any global constant to them still is valid
% we need to fix two DoFs for lambda_0 and s_0 as well
freeDof([1; idxBdSb; idxBdQb]) = false;

%% Record assembling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% direct solve
bigSoln = zeros(NDof,1);

% [~,sys] = memory;
% if sum(freeDof) > 2e6 && sys.PhysicalMemory.Available < 3e10
%     % 2 million DoFs need about 27GB free memory to be solved
%     % 3 million DoFs need about 39GB free memory
%     fprintf('Number of DoF %d exceeds estimated max dofs allowed in memory, stop. \n',sum(freeDof)); 
%     return
% end

bigSoln(freeDof) = bigA(freeDof, freeDof)\Rhs(freeDof);
%% checking error
idxUDoF = (4*NT+3*NF+1:4*NT+3*NF+12*NT)';
uh = bigSoln(idxUDoF);

% q0 starts at (4*NT + 3*NF) + 12*NT + (4*NT+3*NF)= 20*NT+6*NF
idxq0DoF = (20*NT+6*NF+1:20*NT+6*NF+12*NT)';
qh0 = bigSoln(idxq0DoF);

% q has 12*NT+4*NF DoFs
idxqDoF = (20*NT+6*NF+1 : 20*NT+6*NF + 12*NT+6*NF)';
qh = bigSoln(idxqDoF);

idxlambdaDoF = (1:4*NT+3*NF)';
lambdah = bigSoln(idxlambdaDoF);

idxsDoF = (4*NT+3*NF+12*NT+1:4*NT+3*NF+12*NT+4*NT+3*NF)';
sh = bigSoln(idxsDoF);
%%
errL2U = zeros(NT,1);
errL2q = zeros(NT,1);

for p = 1:nQuadBulk
    pxyz = lambdaBulk(p,1)*node(elem(:,1),:) ...
        + lambdaBulk(p,2)*node(elem(:,2),:) ...
        + lambdaBulk(p,3)*node(elem(:,3),:) ...
        + lambdaBulk(p,4)*node(elem(:,4),:);
    
    up = pde.exactu(pxyz);
    
    uhpElem = zeros(NT,3);
    qh0pElem = zeros(NT,3);
    for j = 1:12
        je = j; if j > 6; je = je-6; end
        j1 = edge_to_vert(je,1); j2 = edge_to_vert(je,2);
        
        if j <= 6
            phi_j = lambdaBulk(p,j1)*Dphi(:,:,j2)-lambdaBulk(p,j2)*Dphi(:,:,j1);
        else
            phi_j = lambdaBulk(p,j1)*Dphi(:,:,j2)+lambdaBulk(p,j2)*Dphi(:,:,j1);
        end
        uhpElem = uhpElem + repmat(uh(elem2Vhdof(:,j)),1,3).*phi_j;
        qh0pElem = qh0pElem + repmat(qh(elem2Vhdof(:,j)),1,3).*phi_j;
    end
    
    EpsUminusUh = sum(bsxfun(@times, Eps2elem, up - uhpElem), 2);
    errL2U = errL2U + wBulk(p)*dot(squeeze(EpsUminusUh), up - uhpElem, 2);
    % errL2U = errL2U + wBulk(p)*sum((up - uhElem).^2, 2);
    errL2q = errL2q + wBulk(p)*sum(qh0pElem.^2, 2);
end


errL2U = sqrt(errL2U.*volume);
errL2q = sqrt(errL2q.*volume);

errL2UTotal = norm(errL2U);
errL2QTotal = norm(errL2q);

errorStabTotal = sqrt(qh'*Stab_Wh*qh+lambdah'*Stab_Mh*lambdah);
errorsTotal = sqrt(sh'*Stab_Sh*sh);

%% for visualization
uhElem = zeros(NT,3);
for j = 1:12
    je = j; if j > 6; je = je-6; end
    j1 = edge_to_vert(je,1); j2 = edge_to_vert(je,2);
    
    if j <= 6
        phi_j = 0.25*(Dphi(:,:,j2)-Dphi(:,:,j1));
    else
        phi_j = 0.25*(Dphi(:,:,j2)+Dphi(:,:,j1));
    end
    uhElem = uhElem + repmat(uh(elem2Vhdof(:,j)),1,3).*phi_j;
end
%% Output information
if nargout == 1
    soln = bigSoln;
else
    soln = struct('u',uh,...
                  'uh2elem', uhElem,...
                  'q0',qh0);
    eqn = struct('A',bigA,'b',Rhs,...
                'face',face,...
                'freeDof',freeDof,...
                'bdDofScalar',idxBdFace,...
                'bdDofVec',idxBdFaceVec,...
                'errorU',errL2UTotal, ...
                'errorQ', errL2QTotal,...
                'errorqlambda', errorStabTotal,...
                'errorS', errorsTotal,...
                'errorUelem',errL2U,...
                'NFace',NF);
    info.assembleTime = assembleTime;
end


%% end of function
end