function [bc,mapper] = boundaryEdge3D(p,e,bcind)

%% Usage:
%
% INPUTS:
% mesh --- a struct data contains very rich mesh information.;
% p --- (x,y) coordinates of each node.
% bcind --- [bc1,bc2,bc3,bc4,bc5,bc6] denote the boundary condition on the
%            boundary of a cubic domain, where 1,2,3,4,5,6 denotes the
%            bc1: x = xmin  bc2: x = xmax
%            bc3: y = ymin  bc4: y = ymax
%            bc5: z = zmin  bc6: z = zmax
%            if bc = 1: Dirichlet boundary condition
%               bc = 2: Neumann boundary condition
%               bc = 3: Robin boundary condition
%
% OUTPUTS:
% bc --- indices of global degrees of freedom on the boundary
% mapper --- a vector used to extract out the unknowns from the global
%            degrees of freedom
%
% Remark: this code only handles Dirichlet boundary.
%
% Last Modified: 05/26/2016 by Xu Zhang
% Last Modified: 07/05/2020 by Xu Zhang
%% Form bc
%p = mesh.p;
xmin = min(p(:,1)); xmax = max(p(:,1));
ymin = min(p(:,2)); ymax = max(p(:,2));
zmin = min(p(:,3)); zmax = max(p(:,3));
domain = [xmin,xmax,ymin,ymax,zmin,zmax];

tbc = find(bcind ~= 2);
domain2 = domain(tbc);

p2 = (p(e(:,1),:) + p(e(:,2),:))/2;
p3 = [p2(:,1),p2(:,1),p2(:,2),p2(:,2),p2(:,3),p2(:,3)];
p3 = p3(:,tbc);
np = size(p2,1); mapper = (1:np)';

if isempty(tbc)
    bc = [];    
else
    k = 0;
    for i = 1:length(tbc)
        idtmp = find(abs(p3(:,i) - domain2(i)) < 10*eps);
        nk = length(idtmp);
        id(k+1:k+nk) = idtmp;
        k = k + nk;
    end
    bc = unique(id);
    mapper(bc) = [];
end
