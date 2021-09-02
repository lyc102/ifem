function eqn = mixedPoisson(node,elem,pde,bdflag,option,m)
%% the RT0 - P0 finite element method for mixed Poisson problem
% The weak problem is: Find H_h \in V_h = RT0 and \phi_h \in W_h = P0 such
% that
% (H_h,\psi_h) + (\phi_h ,\div\psi_h) = 0   \forall \psi_h \in V_h
%
% (\div H_h,r_h) = -(div He,r_h) - (\div m_h,r_h) \forall r_h \in W_h

%% Construct Data Structure
[elem2face,face] = dof3face(elem);
locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3]; % ascend ordering
NP = size(node,1);   NT = size(elem,1);  NF = size(face,1);

[Dlambda,volume,elemSign] = gradbasis3(node,elem);

%% Assemble matrix 
Ndof = NF + NT;

%% M.  Mass matrix for RT0 element
% Given a face f_l formed by vertices [i,j,k] the basis for the RT0 element reads as:
% \phi_l = 2(\lambda_i Dlambda_j \times D\lambda_k + \lambda_j Dlambda_k
% \times Dlambda_i + \lambda_k Dlambda_i \times Dlambda_j)
% the corresponding degree of freedom 
% l_f(v) = \int_{f_l} v\cdot n\dd s
ii = zeros(16*NT,1); jj = zeros(16*NT,1); sM = zeros(16*NT,1);
index = 0;
for i = 1:4
    for j = 1:4
        i1 = locFace(i,1); i2 = locFace(i,2); i3 = locFace(i,3);
        j1 = locFace(j,1); j2 = locFace(j,2); j3 = locFace(j,3);
        ii(index+1:index+NT) = double(elem2face(:,i));
        jj(index+1:index+NT) = double(elem2face(:,j));
        % the inner product (\phi_i,\phi_j)
        Mij = 1/5.*volume .* (...
            (1+(i1 == j1)).*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2),...
                                    mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)+...
            (1+(i1 == j2)).*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2),...
                                    mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)+...
            (1+(i1 == j3)).*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2),...
                                    mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)+...
            (1+(i2 == j1)).*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2),...
                                    mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)+...
            (1+(i2 == j2)).*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2),...
                                    mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)+...
            (1+(i2 == j3)).*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2),...
                                    mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)+...
            (1+(i3 == j1)).*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2),...
                                    mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)+...
            (1+(i3 == j2)).*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2),...
                                    mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)+...
            (1+(i3 == j3)).*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2),...
                                    mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2));
        sM(index+1:index+NT) = Mij;
        index = index + NT;
    end
end
M = sparse(ii,jj,sM,NF,NF);
clear ii jj sM Mij;
%M1 = getmassmatvec3(elem2face,volume,Dlambda,'RT0',K);

% B. the transpose of the divergence operator NF \times NT
ii = zeros(4*NT,1); jj = zeros(4*NT,1); sB = zeros(4*NT,1);
index = 0;
% the divergence of the basis \phi_l reads as
% \div\phi_l = 6 Dlambda_i \cdot (Dlambda_j \times Dlambda_k)
for i = 1:4
    i1 = locFace(i,1); i2 = locFace(i,2); i3 = locFace(i,3);
    ii(index+1:index+NT) = double(elem2face(:,i));
    jj(index+1:index+NT) = 1:NT;
    % (div\phi_i, 1) 
    sB(index+1:index+NT) = 6  .* volume .* dot(Dlambda(:,:,i1),mycross(Dlambda(:,:,i2),...
        Dlambda(:,:,i3),2),2);
    index = index + NT;
end
B = sparse(ii,jj,sB,NF,NT);
clear ii jj sB;
%BB = icdmat(double(elem2face),[1 -1 1 -1]);

% C. zero matrix.
C = eye(NT,NT);

A = [M B;B' 10^(-10).*C];
clear C;


%% compute the right hand side
f = zeros(NT,1);
[lambda,weight] = quadpts3(option.fquadorder);
nQuad = size(lambda,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    divHep = pde.f(pxyz);
    f = f - weight(p) .* volume .* divHep;
end
f = f - B' * m;   % B*m = (\div m,r_h)
F = [zeros(NF,1);f];

%% deal with the boundary condition
%[bdNode,bdFace,isBdNode] = findboundary3(elem);
% find boundary face
isDirichlet(elem2face(bdflag(:) == 1)) = true;
isbdFace(isDirichlet) = true;
freeFace = ~isbdFace;
freeDof = find(freeFace == 1);
Dofphi = 1:NT;
Dofphi = Dofphi+NF;
bdFacenum = find(isbdFace == 1);
freeDof = [freeDof,Dofphi];

% deal with the zero boundary condition 
A(bdFacenum,:) = 0;
A(:,bdFacenum) = 0;
F(bdFacenum) = 0;
for i = 1:length(bdFacenum)
    A(bdFacenum(i),bdFacenum(i)) = 1;
end

%% solve the linear equation
ubig = zeros(Ndof,1);
ubig(freeDof) = A(freeDof,freeDof)\F(freeDof);
Hh = ubig(1:NF);

HI = faceinterpolate3(pde.H,node,elem);
%HIvec = HI(elem2face);

errI  = getL2error3RT0(node,elem,pde.H,HI);

errh = getL2error3RT0(node,elem,pde.H,Hh);

L2normH = getL2error3RT0(node,elem,pde.H,zeros(NF,1));

errI = errI / L2normH;
errh = errh / L2normH;

%% error estimates

L2errH = zeros(NT,1); L2normH = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ... 
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    Hp = pde.H(pxyz);
    Hhp = zeros(NT,3);
    for i = 1:4
        i1 = locFace(i,1); i2 = locFace(i,2); i3 = locFace(i,3);
        Hhp = Hhp + 2 * repmat(Hh(elem2face(:,i)),1,3) .* ...
               (lambda(p,i1) .* mycross(Dlambda(:,:,i2),Dlambda(:,:,i3))...
            + lambda(p,i2) .* mycross(Dlambda(:,:,i3),Dlambda(:,:,i1))...
            + lambda(p,i3) .* mycross(Dlambda(:,:,i1),Dlambda(:,:,i2)));
    end
    L2errH = L2errH + weight(p) * sum((Hp - Hhp).^2,2);
        
    L2normH = L2normH + weight(p) * sum(Hp.^2,2);
end
L2errH = volume .* L2errH;
% modify the error
L2errH(isnan(L2errH)) = 0;
L2normH = L2normH .*  volume;
L2errH = sqrt(sum(L2errH));
L2normH = sqrt(sum(L2normH));
eqn = [L2errH/L2normH errI];

