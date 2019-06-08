function uI = edgeinterpolate1(u,node,edge)
%% EDGEINTERPOLATE1 interpolate to the linear edge finite element space.
%
% uI = edgeinterpolate1(u,node,edge) interpolates a given function u into
% the second family linear edge finite element spaces. The coefficient is
% given by the line integral int_e u*t ds and int_e 3*u*t(lambdai-lambdai)
% ds. Simpson rule is used to evaluate the line integrals.
%
% The input u could be a funtional handel or an array of length N (linear
% element) or N+NE (quadratic element).
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
%         uI = edgeinterpolate1(pde.exactu,node,edge);
%         HcurlErr(k) = getHcurlerror3NE1(node,elem,pde.curlu,uI);
%         L2Err(k) = getL2error3NE1(node,elem,pde.exactu,uI);
%         N(k) = length(uI);
%     end
%     figure;
%     showrate2(N,HcurlErr,1,'r-+','||u-u_I||_{curl}',...
%               N,L2Err,1,'b-+','||u-u_I||');
%
%
% See also edgeinterpolate, edgeinterpolate2
%
% <a href="matlab:ifem Maxwelldoc">Maxwell doc</a> Section: Dirichlet
% boundary condition
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = size(node,1);	NE = size(edge,1);
if isnumeric(u)  % allow an array input
    switch length(u)
        case N      % linear element
            uQ = u;
            uQ(N+1:N+NE,:) = (u(edge(:,1),:) + u(edge(:,2),:))/2;
        case N+NE   % quadratic element
            uQ = u;
        case 2*NE   % linear edge element
            uI = u;
            return;
    end
else
    mid = (node(edge(:,1),:) + node(edge(:,2),:))/2;
    uQ = u([node; mid]);
end
edgeVec = node(edge(:,2),:)-node(edge(:,1),:);
uI(1:NE,1) = dot(edgeVec,(uQ(edge(:,1),:)+uQ(edge(:,2),:)+4*uQ(N+1:N+NE,:))/6,2);
uI(NE+1:2*NE,1) = dot(edgeVec,0.5*(uQ(edge(:,1),:)-uQ(edge(:,2),:)),2);