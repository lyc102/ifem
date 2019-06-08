function [sigma,u,eqn,info] = elasticity3mfemP1P0(node,elem,pde,bdFlag,option)


if ~exist('option','var'), option = []; end

time = cputime;  % record assembling time

d = 3; % three dimension
N = size(node,1);
NT = size(elem,1);
Nsigma = 6*N;
Nu = 3*NT;
Ndof = Nsigma + Nu;

%% Assemble matrix
[Dphi,volume] = gradbasis3(node,elem);

%% Mass matrix of linear (P1) element
M = sparse(N,N);
for i = 1:4
    for j = i:4
        ii = double(elem(:,i));
        jj = double(elem(:,j));
        if (j==i)
            M = M + sparse(ii,jj,volume/10,N,N);
        else
            M = M + sparse([ii;jj],[jj;ii],[volume/20; volume/20],N,N);                               
        end                    
    end
end

%% Compliance tensor
lambda = pde.lambda;
mu = pde.mu;
C1 = 1/(2*mu);
C2 = C1*lambda/(2*mu + d*lambda);
% E = mu*(3*lambda+2*mu)/(lambda+mu);
% nu = lambda/(2*(lambda+mu));
A = sparse(C1*eye(6,6) - C2*([1 1 1 0 0 0]'*[1 1 1 0 0 0]));

%% Matrix for (Asigma,tau)
Am = kron(A,M);

%% Div operator
elem2dofsigma = zeros(NT,4,6);
elem2dofsigma(:,:,1) = double(elem);
for k = 2:6
    elem2dofsigma(:,:,k) = elem2dofsigma(:,:,k-1)+N;
end
elemIdx = (1:NT)';
elem2dofu = [elemIdx elemIdx+NT elemIdx+2*NT];
Dx = squeeze(Dphi(:,1,:)).*repmat(volume,1,4);
Dy = squeeze(Dphi(:,2,:)).*repmat(volume,1,4);
Dz = squeeze(Dphi(:,3,:)).*repmat(volume,1,4);
clear Dphi
% u1: dx sigma(1) + dy sigma(4) + dz sigma(6)
% u2: dx sigma(4) + dy sigma(2) + dz sigma(5)
% u3: dx sigma(6) + dy sigma(5) + dz sigma(3)
B = sparse(repmat(elem2dofu(:,1),1,4),elem2dofsigma(:,:,1),Dx,Nu,Nsigma) ...
  + sparse(repmat(elem2dofu(:,1),1,4),elem2dofsigma(:,:,4),Dy,Nu,Nsigma) ...  
  + sparse(repmat(elem2dofu(:,1),1,4),elem2dofsigma(:,:,6),Dz,Nu,Nsigma) ...  
  + sparse(repmat(elem2dofu(:,2),1,4),elem2dofsigma(:,:,4),Dx,Nu,Nsigma) ...
  + sparse(repmat(elem2dofu(:,2),1,4),elem2dofsigma(:,:,2),Dy,Nu,Nsigma) ...  
  + sparse(repmat(elem2dofu(:,2),1,4),elem2dofsigma(:,:,5),Dz,Nu,Nsigma) ...  
  + sparse(repmat(elem2dofu(:,3),1,4),elem2dofsigma(:,:,6),Dx,Nu,Nsigma) ...
  + sparse(repmat(elem2dofu(:,3),1,4),elem2dofsigma(:,:,5),Dy,Nu,Nsigma) ...  
  + sparse(repmat(elem2dofu(:,3),1,4),elem2dofsigma(:,:,3),Dz,Nu,Nsigma);

%% Stabilization
T = auxstructure3(elem);
face2elem = T.face2elem;
elem2face = T.elem2face;
[normal,area,unitNormal] = facenormal(node,T.face);
clear T;
h = sqrt(area);
harea = h.*area;
% elementwise part: [u_in_k']:[u_jn_k']=n_k(i)n_k(j);
% Use Dphi as a scaled outwards normal vector to compute - int_F h[u n'][v n']
C = sparse(Nu,Nu);
for i = 1:3
    for j = i:3
        ii = elem2dofu(:,i);
        jj = elem2dofu(:,j);
        Cij = zeros(NT,1);
        for k = 1:4 % sum of all four faces
            fk = elem2face(:,k);
            Cij = Cij + harea(fk).*((i==j) + unitNormal(fk,i).*unitNormal(fk,j));
        end
        if (j==i)            
            C = C + sparse(ii,jj,-Cij,Nu,Nu);
        else
            C = C + sparse([ii;jj],[jj;ii],repmat(-Cij,1,2),Nu,Nu);                               
        end                    
    end
end
% cross face part
intFaceIdx = (face2elem(:,1) ~= face2elem(:,2));
for i = 1:3
    for j = 1:3
        ii = elem2dofu(face2elem(intFaceIdx,1),i);
        jj = elem2dofu(face2elem(intFaceIdx,2),j);
        Cij = harea(intFaceIdx).*((i==j) + unitNormal(intFaceIdx,i).*unitNormal(intFaceIdx,j));
        C = C + sparse([ii;jj],[jj;ii],repmat(-Cij,1,2),Nu,Nu);                               
    end
end

%% Assemble the right hand side
F = zeros(Ndof,1);
fu = zeros(NT,3);
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isempty(pde.f)
	[lambda,weight] = quadpts3(option.fquadorder);
	nQuad = size(lambda,1);
    for p = 1:nQuad
		% quadrature points in the x-y-z coordinate
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:) ...
             + lambda(p,4)*node(elem(:,4),:);
		fp = pde.f(pxyz);
		fu = fu - weight(p)*fp;
    end
    fu = fu.*repmat(volume,1,3);
end
clear fp
F((Nsigma+1):Ndof,1) = fu(:);

%% Boundary Conditions
if ~exist('bdFlag','var'), bdFlag = []; end
eqn = struct('Am',Am,'B',B,'C',C,'f',F(1:Nsigma),'g',F(Nsigma+1:end));
assembleTime = cputime - time;
%% Solver
bigA = [Am B'; B C];
bigu = bigA\F;
sigma = bigu(1:Nsigma);
u = bigu(Nsigma+1:end);

%% Output information
info.assembleTime = assembleTime; 