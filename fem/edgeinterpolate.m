 function uI = edgeinterpolate(u,node,edge,quadOrder)
%% EDGEINTERPOLATE interpolate to the lowest order edge finite element space.
%
% uI = edgeinterpolate(u,node,edge) interpolates a given function
% u into the lowesr order edge finite element spaces. The coefficient
% is given by the line integral int_e u*t ds. Simpson rule is used to
% evaluate this line integral.
%
% The input u could be a funtional handle of a vector field or an array of length
%  - N (linear element): the grad of this function will be interpolated to ND_0
%  - N+NE (quadratic element): the grad of this function will be interpolated to ND_0
%
% Example
% 
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     maxIt = 4;
%     HcurlErr = zeros(maxIt,1);
%     L2Err = zeros(maxIt,1);
%     N = zeros(maxIt,1);
%     for k = 1:maxIt
%         [node,elem] = uniformbisect3(node,elem);
%         [elem2dof,edge] = dof3edge(elem);
%         pde = Maxwelldata2;
%         uI = edgeinterpolate(pde.exactu,node,edge);
%         HcurlErr(k) = getHcurlerror3ND(node,elem,pde.curlu,uI);
%         L2Err(k) = getL2error3ND(node,elem,pde.exactu,uI);
%         N(k) = length(uI);
%     end
%     figure;
%     showrate2(N,HcurlErr,1,'r-+','||u-u_I||_{curl}',...
%               N,L2Err,1,'b-+','||u-u_I||');
%
% See also edgeinterpolate1, edgeinterpolate2
%
% <a href="matlab:ifem Maxwelldoc">Maxwell doc</a> Section: Dirichlet
% boundary condition
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = size(node,1);	d = size(node,2);  NE = size(edge,1);
edgeVec = node(edge(:,2),:)-node(edge(:,1),:);
if isnumeric(u)  % allow an array input
    switch size(u,1)
        case {N, N+NE}  % continuous linear element or quadratic's gradient
            uI = u(edge(:,2)) - u(edge(:,1));
            % edge bubble on e has no contribution to  int_e grad(u)*t ds
        case NE     % edge element
            uI = u;
        case 1      % u is constant
            uQ = u*ones(N+NE,d);
            uI = dot(edgeVec,(uQ(edge(:,1),:)+uQ(edge(:,2),:)+4*uQ(N+1:N+NE,:))/6,2);
    end   
else % u is a function
    if ~exist('quadOrder','var'), quadOrder = 3; end
    [lambda,weight] = quadpts1(quadOrder);
    nQuad = size(lambda,1);
    K = zeros(NE,d);
    for p = 1:nQuad
        pxyz = lambda(p,1)*node(edge(:,1),:) ...
            + lambda(p,2)*node(edge(:,2),:);
        K = K + weight(p)*u(pxyz);
    end
    uI = dot(edgeVec,K,2);
end