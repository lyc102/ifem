%% SQUAREEDDYCURRENT solves the eddy current equation in a square  using lowest order element.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%% Problem
pde = eddycurrentdata1;
[node,elem] = squaremesh([0,1,0,1],1/16);
bdFlag = setboundary(node,elem,'Dirichlet');

%% Parameters 
maxIt = 4; 
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2 = zeros(maxIt,1); 
errHcurl = zeros(maxIt,1); 

%% Finite Element Method        
for k = 1:maxIt
    % refine mesh    
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    % solve the equation
    [u,edge,eqn] = eddycurrent(node,elem,bdFlag,pde);
    % compute error
    % In 2-D, ND is rotation of RT
    errHcurl(k) = getHdiverrorRT0(node,elem,0,u); 
    errL2(k) = getL2errorNE(node,elem,pde.exactu,u);
    % record information
    N(k) = length(u);
    h(k) = 1./(sqrt(size(node,1))-1);    
end

%% Plot
figure;
showrateh2(h,errHcurl,1,'-*','|| curl (u-u_h)||',...
           h,errL2,1,'k-+','|| u-u_h||');
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|| curl (u-u_h)||'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errHcurl,'%0.5e');