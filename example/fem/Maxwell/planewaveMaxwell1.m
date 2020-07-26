%% PLANEWAVEMAXWELL1 plane wave solutions to Maxwell equations in a cube.
%
% Test Maxwell function can solve problems with complex solution and
% indefinite case.
%
% See also  planewaveMaxwell, planewaveMaxwell2
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all;

%% Generate an initial mesh 
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
bdFlag = setboundary3(node,elem,'Neumann');

%% Parameters
maxIt = 4; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);
energyErrImag = zeros(maxIt,1);
L2ErrImag = zeros(maxIt,1);

%% Get the data of the pde
global d P omega
d = [1, pi/2, pi/2];  % in sphereical coordinate
P = [1, 0, 0];
omega = 1;
r = d(1); theta = d(2); phi = d(3);
d = [r*sin(phi)*cos(theta), r*sin(phi)*sin(theta), r*cos(phi)];
% pde = planewavedataC; % plane wave with complex coefficients
pde = planewavedata; % plane wave with real coefficients and complex solution

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    % solve the equation
%     [u,edge,A,M] = Maxwell(node,elem,HB,pde,bdFlag);
    [u,edge,eqn] =  Maxwell1(node,elem,bdFlag,pde); 
    % compute the error
    uI = edgeinterpolate1(pde.exactu,node,edge);
    L2Err(k) = sqrt(abs(real(u-uI)'*eqn.M*real(u-uI)));    
    energyErr(k) = sqrt(abs(real(u-uI)'*eqn.A*real(u-uI)) + L2Err(k)^2);
    L2ErrImag(k) = sqrt(abs(imag(u-uI)'*eqn.M*imag(u-uI)));
    energyErrImag(k) = sqrt(abs(imag(u-uI)'*eqn.A*imag(u-uI)) + L2ErrImag(k)^2);
%     energyErr(k) = getHcurlerror3ND1(node,elem,pde.curlu,u);
%     L2Err(k) = getL2error3ND1(node,elem,pde.exactu,u);
    N(k) = length(u);
    h(k) = 1./(size(node,1)^(1/3)-1);       
end

%% Plot convergence rates
figure(1); clf; 
showrateh2(h,energyErr,1,'r-+','|| Re(u-u_h)||_A',...
           h,L2Err,1,'b-+','|| Re(u-u_h)||');
figure(2); clf; 
showrateh2(h,energyErrImag,1,'r-+','|| Img(u-u_h)||_A',...
           h,L2ErrImag,1,'b-+','|| Img(u-u_h)||');