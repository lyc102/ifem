function eta = estimateresidualWG(node,elem,u,Du,pde)

% Du: graident of CR part

NT = size(elem,1);

if ~exist('Du','var')
   [Du,area] = gradu(node,elem,u); 
end

%% data structure
T = auxstructure(elem);
% neighbor = T.neighbor;
edge2elem = T.edge2elem;
elem2edge = T.elem2edge;
edge = T.edge;
% elem2dof = NT + T.elem2edge;
clear T

%% Weight of diffusion (scalar diffusion and piecewise constant)
if exist('pde','var') && isfield(pde,'d') && ~isempty(pde.d)
    if isreal(pde.d)
        AT = pde.d;            % d is an array
    else                            % d is a function
        center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
        AT = pde.d(center);              
    end
%     Aemax = max(AT(edge2elem,[],2));
    Aemin = min(AT(edge2elem(:,1:2)),[],2);
else
    Aemin = 1;
end      

%% Modification for singular nodes
if exist('pde','var') && isfield(pde,'singularnode') 
   singularNode =  pde.singularnode;
   idx = (elem(:,1) == singularNode) | (elem(:,2) == singularNode) | ... 
         (elem(:,3) == singularNode);     
   ATomegax = AT(idx);
   Amax = max(ATomegax);
   edgeIdx(elem2edge(idx,:)) = true;
   LambdaTbar(edge2elem(edgeIdx,1:2)) = 1;
%    Amin = min(ATomegax);
%    LambdaT = ATomegax/Amin;
   LambdaTbar(idx) = Amax*ATomegax.^(-1);
%    Lambdae = max(LambdaT(edge2elem(edgeIdx,1:2),[],2));
   LambdaEbar = max(LambdaTbar(edge2elem(edgeIdx,1:2)),[],2);
end

%% Data oscillation
osc = zeros(NT,1);
if isfield(pde,'f') && isnumeric(pde.f) && (pde.f==0)
    pde.f = [];
end
if isfield(pde,'f') && ~isempty(pde.f)
    fc = pde.f(center);
    [lambda,weight] = quadpts(3);
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        osc = osc + weight(p)*(fp-fc).^2;
    end
    osc = osc.*(area.^4)./AT;
end

%% Jump of tangential component of the gradient
edgeVec = node(edge(:,2),:) - node(edge(:,1),:);
t1 = edge2elem(:,1);
t2 = edge2elem(:,2);
jumpatEdge = Aemin.*dot(Du(t1,:) - Du(t2,:),edgeVec,2).^2; % edge length is included
if exist('LambdaEbar','var')
    jumpatEdge(edgeIdx) = LambdaEbar.*jumpatEdge(edgeIdx);
end
% merge to elementwise estimator
edgeJumpb = jumpatEdge(elem2edge(:,1)) ...
          + jumpatEdge(elem2edge(:,2)) ...
          + jumpatEdge(elem2edge(:,3));
edgeJumpo = zeros(NT,1);
% edge jump of Duo.
% geometry quantity
% center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
% mid1 = (node(elem(:,2),:) + node(elem(:,3),:))/2;
% mid2 = (node(elem(:,3),:) + node(elem(:,1),:))/2;
% mid3 = (node(elem(:,1),:) + node(elem(:,2),:))/2;
% d1 = mid1 - center;
% d2 = mid2 - center;
% d3 = mid3 - center;
% ct2 = 3./sum(d1.^2 + d2.^2 + d3.^2,2);
% uc = 2*((u(elem2dof(:,1)) + u(elem2dof(:,2)) + u(elem2dof(:,3)))/3 - u(1:NT));
% uc is a high order term and thus can be skipped
% j1 = uc.*ct2.*dot(d1,ve1,2);
% j2 = uc.*ct2.*dot(d2,ve2,2);
% j3 = uc.*ct2.*dot(d3,ve3,2);
% % the sign is already in ve which is different in neighbor element
% edgeJumpo = (j1 + j1(neighbor(:,1))).^2 ...
%           + (j2 + j2(neighbor(:,2))).^2 ...
%           + (j3 + j3(neighbor(:,3))).^2;

%% Jump of normal component of the flux
% For piecewise constant diffusion coefficient, the normal jump of flux can
% be skipped.
%
% ve1 = node(elem(:,3),:)-node(elem(:,2),:);
% ve2 = node(elem(:,1),:)-node(elem(:,3),:);
% ve3 = node(elem(:,2),:)-node(elem(:,1),:);
% area = 0.5*abs(-ve3(:,1).*ve2(:,2)+ve3(:,2).*ve2(:,1));
% % edge jump of Dub
% edgeJumpb = dot((Du-Du(neighbor(:,1),:)),ve1,2).^2 ...
%           + dot((Du-Du(neighbor(:,2),:)),ve2,2).^2 ...
%           + dot((Du-Du(neighbor(:,3),:)),ve3,2).^2;      


%% Residual type error estimator
eta = (osc + 0.5*edgeJumpb + 0.5*edgeJumpo).^(1/2);