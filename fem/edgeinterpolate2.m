function uI = edgeinterpolate2(exactu,node,edge,face,face2edge)
%% EDGEINTERPOLATE2 interpolate to the quadratic (1st type) edge finite element space.
%
% uI = edgeinterpolate(u,node,edge) interpolates a given function
% u into the lowesr order edge finite element spaces. The coefficient
% is given by the line integral int_e u*t ds. Simpson rule is used to
% evaluate this line integral.
%
% The input u could be a funtional handel or an array of length N (linear
% element) or N+NE (quadratic element).
%
% Example
% 
%   node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%   elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%   maxIt = 3;
%   HcurlErr = zeros(maxIt,1);
%   L2Err = zeros(maxIt,1);
%   N = zeros(maxIt,1);
%   for k = 1:maxIt
%       [node,elem] = uniformbisect3(node,elem);
%       [elem2dof,T] = dof3NE2(elem);
%       pde = Maxwelldata2;
%       uI = edgeinterpolate2(pde.exactu,node,T.edge,T.face,T.face2edge);
%       HcurlErr(k) = getHcurlerror3NE2(node,elem,pde.curlu,uI);
%       L2Err(k) = getL2error3NE2(node,elem,pde.exactu,uI);
%       N(k) = length(uI);
%   end
%   figure(1)
%   r1 = showrate(N,HcurlErr,1,'r-+');
%   hold on
%   r2 = showrate(N,L2Err,1,'b-+');
%   legend('||u-u_I||_{curl}',['N^{' num2str(r1) '}'],...
%          '||u-u_I||',['N^{' num2str(r2) '}'],'LOCATION','Best');
%
% See also edgeinterpolate, edgeinterpolate1
%
% <a href="matlab:ifem Maxwelldoc">Maxwell doc</a> Section: Dirichlet
% boundary condition
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Evaluate function at vertices and middle points
if isnumeric(exactu)
    uQ = exactu;
else
    mid = (node(edge(:,1),:) + node(edge(:,2),:))/2;
    uQ = exactu([node; mid]);
end

%% Edge dof coefficients
N = size(node,1);	NE = size(edge,1);    NF = size(face,1);
edgeVec = node(edge(:,2),:)-node(edge(:,1),:);
uI = zeros(2*(NE+NF),1);
uI(1:NE,1) = dot(edgeVec,(uQ(edge(:,1),:)+uQ(edge(:,2),:)+4*uQ(N+1:N+NE,:))/6,2);
uI(NE+1:2*NE,1) = dot(edgeVec,0.5*(uQ(edge(:,1),:)-uQ(edge(:,2),:)),2);

%% Face dof coefficients
eik = node(face(:,3),:)-node(face(:,1),:);
eij = node(face(:,2),:)-node(face(:,1),:);
uquadpts = uQ(N+face2edge(:,1),:)+uQ(N+face2edge(:,2),:)+uQ(N+face2edge(:,3),:);
lf = [4*dot(eik,uquadpts,2) 4*dot(eij,uquadpts,2)];
face2edgeDofValue = [uI(face2edge(:,1:3)) uI(face2edge(:,1:3)+NE)];
localMatrix = [4 8 4 -4  0 4; 8 4 -4 0 -4 4]';
lf = lf - face2edgeDofValue*localMatrix;
uI((2*NE+1):end) = [2*lf(:,1) - lf(:,2); 2*lf(:,2) - lf(:,1)]/3;
