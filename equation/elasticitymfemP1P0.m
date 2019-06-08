function [sigma,u,eqn,info] = elasticitymfemP1P0(node,elem,pde,bdFlag,option)


if ~exist('option','var'), option = []; end

time = cputime;  % record assembling time

d = 2; % two dimensions
N = size(node,1);
NT = size(elem,1);
Nsigma = 3*N;
Nu = 2*NT;
Ndof = Nsigma + Nu;

%% Assemble matrix
[Dphi,area] = gradbasis(node,elem);

%% Mass matrix of linear (P1) element
M = sparse(N,N);
for i = 1:3
    for j = i:3
        ii = double(elem(:,i));
        jj = double(elem(:,j));
        if (j==i)
            M = M + sparse(ii,jj,area/6,N,N);
        else
            M = M + sparse([ii;jj],[jj;ii],[area/12; area/12],N,N);                               
        end                    
    end
end

%% Compliance tensor
lambda = pde.lambda;
mu = pde.mu;
C1 = 1/(2*mu);
C2 = lambda/(2*mu + d*lambda);
% E = mu*(3*lambda+2*mu)/(lambda+mu);
% nu = lambda/(2*(lambda+mu));
A = sparse(C1*(eye(3,3) - C2*([1 1 0]'*[1 1 0])));

%% Matrix for (Asigma,tau)
Am = kron(A,M);

%% Div operator
elem2dofsigma = zeros(NT,3,3);
elem2dofsigma(:,:,1) = double(elem);
for k = 2:3
    elem2dofsigma(:,:,k) = elem2dofsigma(:,:,k-1)+N;
end
elemIdx = (1:NT)';
elem2dofu = [elemIdx elemIdx+NT];
Dx = squeeze(Dphi(:,1,:)).*repmat(area,1,3);
Dy = squeeze(Dphi(:,2,:)).*repmat(area,1,3);
clear Dphi
% u1: dx sigma(1) + dy sigma(3)
% u2: dx sigma(3) + dy sigma(2)
B = sparse(repmat(elem2dofu(:,1),1,3),elem2dofsigma(:,:,1),Dx,Nu,Nsigma) ...
  + sparse(repmat(elem2dofu(:,1),1,3),elem2dofsigma(:,:,3),Dy,Nu,Nsigma) ...  
  + sparse(repmat(elem2dofu(:,2),1,3),elem2dofsigma(:,:,3),Dx,Nu,Nsigma) ...
  + sparse(repmat(elem2dofu(:,2),1,3),elem2dofsigma(:,:,2),Dy,Nu,Nsigma);

%% Stabilization
T = auxstructure(elem);
edge2elem = T.edge2elem;
elem2edge = T.elem2edge;
[normal,edgeLength,unitNormal] = edgenormal(node,T.edge);
clear T;
harea = edgeLength.^2;
% elementwise part: [u_in_k']:[u_jn_k']=n_k(i)n_k(j);
% Use Dphi as a scaled outwards normal vector to compute - int_F h[u n'][v n']
C = sparse(Nu,Nu);
for i = 1:2
    for j = i:2
        ii = elem2dofu(:,i);
        jj = elem2dofu(:,j);
        Cij = zeros(NT,1);
        for k = 1:3 % sum of all four faces
            fk = elem2edge(:,k);
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
for i = 1:2
    for j = 1:2
        ii = elem2dofu(edge2elem(:,1),i);
        jj = elem2dofu(edge2elem(:,2),j);
        Cij = harea.*((i==j) + unitNormal(:,i).*unitNormal(:,j));
        C = C + sparse([ii;jj],[jj;ii],repmat(-Cij,1,2),Nu,Nu);                               
    end
end

%% Assemble the right hand side
F = zeros(Ndof,1);
fu = zeros(NT,2);
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
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
		fu = fu - weight(p)*fp;
    end
    fu = fu.*repmat(area,1,2);
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