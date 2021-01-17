function [soln,eqn,info] = DivCurl3PDWG(node,elem,bdFlag,pde,option,varargin)
%% DivCurl3PDWG div-curl system equation: 
%        primal-dual weak-Galerkin in 3-D to solve the following system.
%       div(epsilon*u)=f  in \Omega, 
%       curl(u) = g       in \Omega, 
%       Dirichlet boundary condition epsilon*u \cdot n =g_n on \Gamma_0, 
%
% The code is vectorized with no loops. Data structure follows iFEM tradition.
%
% Reference: A New Numerical Method for Div-Curl Systems with Low Regularity Assumptions
% S. Cao, C. Wang, J. Wang, https://arxiv.org/abs/2101.03466
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end
if ~isfield(option,'inverseh'), option.inverseh = true; end
if ~isfield(option,'dquadorder'), option.dquadorder = 3; end

[lambda,w] = quadpts3(option.dquadorder); % quadrature order is 1 or 2
nQuad = size(lambda,1);

if ~exist('bdFlag','var')
    bdFlag = setboundary3(node,elem,'Dirichlet');
end


%% Mesh related geometric quantities

% local edge to face indices, first two edges' local indices on a face
loc_e2f = [4 5; 2 3; 1 3; 1 2];
% local face to vertices indices, positive oriented
loc_f2v = [2 3 4; 1 4 3; 1 2 4; 1 3 2];

loc_e2v = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

[elem2face,face] = dof3face(elem);
NT = size(elem,1);
NF = size(face,1); 
elem2intdof = reshape(1:3*NT, 3, NT)';

[Dphi,volume] = gradbasis3(node,elem);

hK = (6*volume).^(1/3);

elem2ve = zeros(NT,3,6);
for e = 1:6 % six edges
   elem2ve(:,:,e) = node(elem(:,loc_e2v(e,2)),:)-node(elem(:,loc_e2v(e,1)),:);
end

elem2ve = elem2ve./repmat(sqrt(sum(elem2ve.^2,2)),1,3); % normalize for basis

faceNormal = cross(node(face(:,2),:) - node(face(:,1),:),...
    node(face(:,3),:) - node(face(:,2),:),2);
faceArea = 0.5*sqrt(sum(faceNormal.^2,2));
faceArea2elem = faceArea(elem2face);
faceArea2elemTotal = sum(faceArea2elem,2);

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

%% assembling
% DoF order: # of Dofs
% \lambda_h: NT + NF, local 1+4
% u_h: 3*NT
% s_h: s_0 + s_b;  NT + NF (simply-connected case), local 1+4
% q_h: q_0 + q_b; 3*NT + 2*NF (for H(curl) vector field)
% assembling indices following locally within each group

%% equation 1 for u_h and lambda_h

blockUlambda = sparse(4*NT+NF,4*NT+NF);  %#ok<NASGU>
% (u_h, \epsilon\nabla_w phi)
% + \sum h_T^{-1} <\lambda_0 - \lambda_b, \varphi_0 - \varphi_b>_{boundary T}

%% (u_h, \epsilon\nabla_w phi_{face})
% locally the basis for P_0(T)^3 is edge [1, 2], [1, 3], and [1, 4]
% size (NT+NF, 3*NT)
C = sparse(NF,3*NT); %(u_h, \epsilon\nabla_w phi_{face})

for j = 1:3 % constant vector space spanned by 3 edge vectors
    for i = 1:4 % face 1 to 4
        % Cij = (q_j, eps \nabla_w \phi_{F_i})
        % \nabla_w \phi_{F_i} = -3\nabla \phi_i
        
        Cij = dot(elem2ve(:,:,j),-3*EpsDphi(:,:,i),2).*volume;
        
        ii = double(elem2face(:,i));
        jj = elem2intdof(:,j);
        C = C + sparse(ii,jj,Cij,NF,3*NT);
    end
end

Ugradphi = [sparse(NT, 3*NT); C];


%% \sum h_T^{-1} <\lambda_0 - \lambda_b, \varphi_0 - \varphi_b>_{boundary T}
% this is essetially a mass matrix stab
Stab_Mh = sparse(NT+NF,NT+NF); %#ok<NASGU>

L0ij = faceArea2elemTotal./hK;
Lambda00 = sparse(1:NT, 1:NT, L0ij, NT, NT); %% <\lambda_0, \lambda_0>

Lambda0b = sparse(NT, NF); %% <\lambda_b, \lambda_0>
Lambdabb = sparse(NF, NF); %% <\lambda_b, \lambda_b>
for j = 1:4 % face 1 to 4
    jj = double(elem2face(:,j));
    Lbj = faceArea2elem(:,j)./hK; % = <\lambda_{F_j}, \lambda_0>
    Lambda0b = Lambda0b + sparse(1:NT, jj, Lbj, NT, NF);
    Lambdabb = Lambdabb + sparse(jj, jj, Lbj, NF, NF);
end

Stab_Mh = [Lambda00 -Lambda0b; -Lambda0b' Lambdabb]; 
% Stablization on M_h should have null space dimension 1 (global constant)

%% \sum h_T <s_0 - s_b, r_0 - r_b>_{boundary T}

Stab_Sh = sparse(NT+NF,NT+NF); %#ok<NASGU>

S0ij = faceArea2elemTotal.*hK;
S00 = sparse(1:NT, 1:NT, S0ij, NT, NT); %% <s_0, s_0>


Sb0 = sparse(NT, NF); %% <s_b, s_0>
Sbb = sparse(NF, NF); %% <s_b, s_b>
for j = 1:4 % face 1 to 4
    jj = double(elem2face(:,j));
    Sbj = faceArea2elem(:,j).*hK; % = <\lambda_{F_j}, \lambda_0>
    Sb0 = Sb0 + sparse(1:NT, jj, Sbj, NT, NF);
    Sbb = Sbb + sparse(jj, jj, Sbj, NF, NF);
end

Stab_Sh = [S00 -Sb0; -Sb0' Sbb]; 
% Stablization on s_h should have null space dimension 1 (global constant)

%% Equation 1 for u_h and q_h = {q_0, q_b}

%%  (u_h, \epsilon\nabla_w \times psi_{face})
% locally the basis for P_0(T)^3 is the edge vector representing 
% the edge [1, 2], [1, 3], and [1, 4] (local indexing of iFEM tradition)
% locally the basis for the tangential component are 
% the edge 1 and 2 vectors for the triangular face (local indexing on this face)
% size: (3*NT+2*NF, 3*NT)

D = sparse(2*NF,3*NT); 
%(u_h, \nabla_w \times psi_{face}) with 2 basis vec at each face

for i = 1:4 % face 1 to 4
    for j = 1:3 % constant vector space spanned by 3 edge vectors
        for e = 1:2 % two basis vectors on each face
            % Dije = (q_j, eps \nabla_w \times \psi_{F_i,e})
            % \nabla_w \times [0, \psi_{F_i}] =
            % (rotation of \psi_{F_i} counter-clockwise by pi/2)*|F_i|/|K|
            % = 3 \psi_{F_i} \times \nabla phi_i
            
            curlpsiFe = 3*mycross(elem2ve(:,:,loc_e2f(i,e)), Dphi(:,:,i));
            Dije = dot(elem2ve(:,:,j),curlpsiFe,2).*volume;
            
            % k-th (global face indexing) face has 2 DoFs in W_h: 2k-1, and 2k
            ii = 2*double(elem2face(:,i)) - e + 1;
            jj = elem2intdof(:,j);
            D = D + sparse(ii,jj, Dije, 2*NF,3*NT);
        end
    end
end

Ucurlpsi = [sparse(3*NT, 3*NT); D];

%%  (\psi, \epsilon\nabla_w s_h)
Psigrads = sparse(3*NT+2*NF, NT+NF);


E = sparse(3*NT,NF);
for i = 1:3 % constant vector space spanned by 3 edge vectors W_{h,0} 
    for j = 1:4 % face 1 to 4
        % Eij = (\psi_{K,i}, eps \nabla_w \phi_{F_j})
        % \nabla_w \phi_{F_j} = -3\nabla \phi_j
        
        Eij = dot(elem2ve(:,:,i),-3*EpsDphi(:,:,j),2).*volume;
        
        ii = elem2intdof(:,i);
        jj = double(elem2face(:,j));
        E = E + sparse(ii,jj,Eij,3*NT,NF);
    end
end

%%% E is C's transpose, the code above is for booking-keeping
E = C';
Psigrads(1:3*NT,NT+1:end) = E;

%% \sum h_T^{-1} <(q_0-q_b)\times n, (\psi_0-\psi_b)\times n>_{boundary T}


% \sum h_T^{-1} <q_0\times n, \psi_0 \times n>_{boundary T}
WS00 = sparse(3*NT, 3*NT);
for i = 1:3 % constant vector space spanned by 3 edge vectors W_{h,0} 
    for j = 1:3 % constant vector space spanned by 3 edge vectors W_{h,0} 
        S00ijface = 0;
        for f = 1:4
            psi0crossn = mycross(elem2ve(:,:,i), Dphi(:,:,f)); % test
            q0crossn = mycross(elem2ve(:,:,j), Dphi(:,:,f)); % trial
            
            S00ijface = S00ijface + 9*dot(q0crossn,psi0crossn,2)...
                .*volume.^2./faceArea2elem(:,f);
        end
        S00ijface = S00ijface./hK;
        ii = elem2intdof(:,i);
        jj = elem2intdof(:,j);
        WS00 = WS00 + sparse(ii,jj,S00ijface,3*NT, 3*NT);
    end
end

% \sum h_T^{-1} <q_b\times n, \psi_0 \times n>_{boundary T}
WSb0 = sparse(3*NT, 2*NF);
for i = 1:3 % constant vector space spanned by 3 edge vectors W_{h,0} 
    for j = 1:4 % 4 face for W_{h,b} each basis has support only on F_j
        for e = 1:2 % each face has 2 DoFs
            psi0crossn = mycross(elem2ve(:,:,i), Dphi(:,:,j)); % test
            qbcrossn = mycross(elem2ve(:,:,loc_e2f(j,e)), Dphi(:,:,j)); % trial
            
            Sb0ijface = 9*dot(qbcrossn,psi0crossn,2)...
                .*volume.^2./faceArea2elem(:,j)./hK;
            
            ii = elem2intdof(:,i);
            % k-th (global face indexing) face has 2 DoFs
            jj = 2*double(elem2face(:,j)) - e + 1;
            WSb0 = WSb0 + sparse(ii,jj,Sb0ijface,3*NT,2*NF);
        end
    end
end

% \sum h_T^{-1} <q_b\times n, \psi_b \times n>_{boundary T}
WSbb = sparse(2*NF, 2*NF);
for f = 1:4 
    % this term is non-zero only when both trial and test
    % are associated with the same face
    for i = 1:2 % each face has 2 DoFs
        for j = 1:2
            psibcrossn = mycross(elem2ve(:,:,loc_e2f(f,i)), Dphi(:,:,f));
            qbcrossn = mycross(elem2ve(:,:,loc_e2f(f,j)), Dphi(:,:,f));
            
            Sbbijface = 9*dot(qbcrossn,psibcrossn,2)...
                .*volume.^2./faceArea2elem(:,f)./hK;
            
            % k-th (global face indexing) face has 2 DoFs
            ii = 2*double(elem2face(:,f)) - i + 1;
            jj = 2*double(elem2face(:,f)) - j + 1;
            WSbb = WSbb + sparse(ii,jj,Sbbijface,2*NF,2*NF);
        end
    end
end

Stab_Wh = [WS00 -WSb0; -WSb0' WSbb]; 
% Stab_Wh should have null space dimension 
% = span{gradient of continuous P_1 Langrange} = # node -1

%% overall stiffness matrix for the system
blockUlambda = [Stab_Mh Ugradphi; Ugradphi' sparse(3*NT,3*NT)];
blockOffDiag = sparse(4*NT+3*NF, 4*NT+NF); 
blockOffDiag(NT+NF+1:end, NT+NF+1:end) = Ucurlpsi;
blockQs = [-Stab_Sh Psigrads'; Psigrads Stab_Wh];
bigA = [blockUlambda blockOffDiag'; blockOffDiag blockQs];

%% right hand sides

%% -(f,\varphi_0)
ft = zeros(NT,1);

for p = 1:nQuad
		% quadrature points in the x-y-z coordinate
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:) ...
             + lambda(p,4)*node(elem(:,4),:);
		fp = pde.f(pxyz);
        ft = ft + w(p)*fp;
end

Fphi0 = -accumarray((1:NT)',ft.*volume,[NT 1]);

%% \sum_{T} <eps u\cdot n, \varphi_b>_{\partial T} near \partial \Omega
GNphib = zeros(NF,1); % <eps u\cdot n, \varphi_F>_{\partial T}
% update: face quadorder increased

if ~isfield(option,'gNquadorder')
    option.gNquadorder = 5;   % default order exact for linear gN
end
[lambdagN, weightgN] = quadpts(option.gNquadorder);               % linear bases
nQuadgN = size(lambdagN,1);

for i = 1:4 % 4 faces
    idxBdElem = find(bdFlag(:,i) > 0); % if an element has a boundary face 

    if ~isempty(idxBdElem)
        bdNormal = -3*repmat(volume(idxBdElem,:),1,3).*Dphi(idxBdElem,:,:); 
        % scaled normal with face area built-in
        loc_face = loc_f2v(i,:);
        
        EpsUdotN = zeros(size(idxBdElem,1),1);
        
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxyz = lambdagN(pp,1)*node(elem(idxBdElem,loc_face(1)),:) ...
                + lambdagN(pp,2)*node(elem(idxBdElem,loc_face(2)),:) ...
                + lambdagN(pp,3)*node(elem(idxBdElem,loc_face(3)),:);
            
            
            u_bd = pde.exactu(ppxyz);
            EpsU = sum(bsxfun(@times, Eps2elem(idxBdElem,:,:), u_bd), 2);
            EpsUdotN = EpsUdotN + ...
                weightgN(pp)*dot(squeeze(EpsU), bdNormal(:,:,i), 2);
               
        end
        
        GNphib = GNphib + accumarray(elem2face(idxBdElem,i), EpsUdotN, [NF 1]);
    end
end

%% (g, \psi_0)
Gpsi0 = zeros(3*NT,1);

for i = 1:3 % 3 vector basis on each element for \psi_0
    gtpsi0 = zeros(NT,1);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:) ...
            + lambda(p,4)*node(elem(:,4),:);
        gp = pde.g(pxyz);
        
        gtpsi0 = gtpsi0 + w(p)*dot(gp, elem2ve(:,:,i), 2);
    end
      
    Gpsi0 = Gpsi0 + accumarray(elem2intdof(:,i), gtpsi0.*volume, [3*NT 1]);
end
%% RHS
Rhs = [Fphi0; GNphib]; % test: phi_0, phi_b
Rhs = [Rhs; zeros(3*NT, 1)]; % test: v
Rhs = [Rhs; zeros(NT, 1)]; % test: r_0
Rhs = [Rhs; zeros(NF, 1)]; % test: r_b
Rhs = [Rhs; Gpsi0]; % test: psi_0
Rhs = [Rhs; zeros(2*NF, 1)]; % test: psi_b

%% set up free dof
% DoF order: # of Dofs
% prefix with "idx": global indices
% the DoF variable itself is boolean (fastest on MATLAB)
% \lambda_h: NT + NF, local 1+4 (no BC)
% u_h: 3*NT (no BC)
% s_h: s_0 + s_b;  NT + NF (simply-connected case), s_b = 0 on Gamma_0
% q_h: q_0 + q_b; 3*NT + 2*NF (for H(curl) vector field) q_b\cross n =0
% assembling indices following locally within each group

isBdFace = false(NF,1);
isBdFace(elem2face(bdFlag(:) > 0)) = true;
idxBdFace = find(isBdFace);
idxBdFaceVec = [2*idxBdFace-1 2*idxBdFace]';
idxBdFaceVec = idxBdFaceVec(:);

freeDof = true(8*NT+4*NF,1);
idxBdSb = idxBdFace + 5*NT+NF; % boundary DoF indices for s_b
idxBdQb = idxBdFaceVec + 8*NT+2*NF; % boundary DoF indices for q_b


%%%%% fixing 1 in lambda_0 is enough, the others are optional %%%%%
%%% current config according to the paper: fixing 1 in lambda_0
%%% bd dofs of s and q
freeDof([1; idxBdSb; idxBdQb]) = false; 

%%%%% fixing bd dofs for q
% freeDof([1; idxBdQb]) = false; 

%%%%% fixing bd dofs for s
% freeDof([1; idxBdSb]) = false;

idxIntBdFace = [];
% old way of handling cavity (second betti number not zero)
if any(elem2face(bdFlag(:) == 2)) && false 
    disp('Modifying cavity dofs')
    isIntBdFace = false(NF,1);
    isIntBdFace(elem2face(bdFlag(:) == 2)) = true;
    idxIntBdFace = find(isIntBdFace);
    idxIntBdSb = idxIntBdFace + 5*NT+NF;  % boundary DoF indices for s_b on cavity
    freeDof(idxIntBdSb) = true;
end

%% Record assembling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% direct solve
NDof = 8*NT+4*NF;

bigSoln = zeros(NDof,1);

[~,sys] = memory;
if sum(freeDof) > 2e6 && sys.PhysicalMemory.Available < 3e10
    % 2 million DoFs need about 27GB free memory to be solved
    % 3 million DoFs need about 39GB free memory
    fprintf('Number of DoF %d exceeds estimated max dofs allowed in memory, stop. \n',sum(freeDof)); 
    return
end
bigSoln(freeDof) = bigA(freeDof, freeDof)\Rhs(freeDof);

%% post-processing for cavity
% TO-DO: current code allows only 1 cavity
if any(elem2face(bdFlag(:) == 2))
    disp('Post-processing cavity dofs')
    intBdFaceVec = false(size(freeDof,1), 1);
    isIntBdFace = false(NF,1);
    isIntBdFace(elem2face(bdFlag(:) == 2)) = true;
    idxIntBdFace = find(isIntBdFace);
    idxIntBdSb = idxIntBdFace + 5*NT+NF;  % boundary DoF indices for s_b on cavity
    intBdFaceVec(idxIntBdSb) = true;
    
    AintBdFaceVec = bigA*intBdFaceVec;   
    sbConst = dot((Rhs - bigA*bigSoln),AintBdFaceVec)/dot(AintBdFaceVec,AintBdFaceVec);
    fprintf('s_b on the interior boundary is %6.4g \n',sbConst);
    bigSoln = bigSoln + sbConst*intBdFaceVec;
end


%% computing error
errL2U = zeros(NT,1);

idxUDoF = (NT+NF+1:NT+NF+3*NT)';
uh = bigSoln(idxUDoF);
uhElem = zeros(NT,3);

idxq0DoF = (5*NT+2*NF+1:5*NT+2*NF+3*NT)';
qh0 = bigSoln(idxq0DoF);
qh0Elem = zeros(NT,3);

idxqbDoF = (8*NT+2*NF+1:8*NT+4*NF)';
qhb = bigSoln(idxqbDoF);

idxqDoF = (5*NT+2*NF+1:8*NT+4*NF)';
qh = bigSoln(idxqDoF);

idxlambdaDoF = (1:NT+NF)';
lambdah = bigSoln(idxlambdaDoF);

lambdah0 = bigSoln((1:NT)');
lambdahb = bigSoln((NT+1:NT+NF)');

idxsDoF = (4*NT+NF+1:5*NT+2*NF)';
sh = bigSoln(idxsDoF);
sh0 = bigSoln((4*NT+NF+1:5*NT+NF)');
shb = bigSoln((5*NT+NF+1:5*NT+2*NF)');

for k = 1:3 % three bases in each K       
    uhElem = uhElem + repmat(uh(elem2intdof(:,k)),1,3).*elem2ve(:,:,k);
    qh0Elem = qh0Elem + repmat(qh0(elem2intdof(:,k)),1,3).*elem2ve(:,:,k);
end


[lambda,w] = quadpts3(option.dquadorder+2); % quadrature order is 1 or 2
nQuad = size(lambda,1);
for p = 1:nQuad
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    
    up = pde.exactu(pxyz);
    EpsUminusUh = sum(bsxfun(@times, Eps2elem, up - uhElem), 2);
    errL2U = errL2U + w(p)*dot(squeeze(EpsUminusUh), up - uhElem, 2);
end


errL2U = sqrt(errL2U.*volume);
errL2UTotal = norm(errL2U);

% idxfreeFace = find(~isBdFace);
% idxFreeFaceVec = [2*idxfreeFace-1 2*idxfreeFace]';
% idxFreeFaceVec = idxFreeFaceVec(:);
% errorStabTotal = sqrt(qh(idxFreeFaceVec)'*...
%     Stab_Wh(idxFreeFaceVec,idxFreeFaceVec)*qh(idxFreeFaceVec)...
%     +lambdah'*Stab_Mh*lambdah);
errorStabTotal = sqrt(qh'*Stab_Wh*qh+lambdah'*Stab_Mh*lambdah);
errorsTotal = sqrt(sh'*Stab_Sh*sh);



%% Output information
if nargout == 1
    soln = bigSoln;
else
    soln = struct('u',uh,...
                  'uh2elem', uhElem,...
                  'q0',qh0,...
                  'qb',qhb,...
                  'lambda0',lambdah0,...
                  'lambdab',lambdahb,...
                  's0',sh0,...
                  'sb',shb);
    eqn = struct('A',bigA,'b',Rhs,...
                'face',face,...
                'freeDof',freeDof,...
                'bdDofScalar',idxBdFace,...
                'bdDofVec',idxBdFaceVec,...
                'bdDofCavity',idxIntBdFace,...
                'errorU',errL2UTotal, ...
                'errorqlambda', errorStabTotal,...
                'errorS', errorsTotal,...
                'errorUelem',errL2U,...
                'NFace',NF);
    info.assembleTime = assembleTime;
end


%% end of function
end