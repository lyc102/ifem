%% FOURCURL3P2TEST
%
% Use the quadraticNedelec elements to approximate the velocity u and
% the stream function w. 
%
% We solve the following equations:
%
%  -w + curl curl u = 0 
%  curl curl w + u  = f
%
% with Dirichlet boundary condition
%
% u\times n = (curl u )\times n  = 0  on \partial \Omega
%
% Please check fourCurl3doc for details.
%
% Lin Zhong, May, 2013. Clean up. Long Chen, Oct 21, 2018
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all; 
clear variables;

%% Problem Setting
pde = fourCurl3data;
[node,elem] = cubemesh([-1,1,-1,1,-1,1],2);
bdFlag = setboundary3(node,elem,'Dirichlet');

%% Parameters
maxIt = 3;
errwIwhL2 = zeros(maxIt,1); 
errwIwhHcurl = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);
errwwhHcurl = zeros(maxIt,1);
erruIuhL2 = zeros(maxIt,1); 
erruIuhHcurl = zeros(maxIt,1);
erruL2 = zeros(maxIt,1);
erruuhHcurl = zeros(maxIt,1);
assembleTime = zeros(maxIt,1);
solverTime = zeros(maxIt,1);
N = zeros(maxIt,1);
h = zeros(maxIt,1);
option = [];

for k = 1:maxIt
    % refine grid        
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    % solve the equation
    [w,u,eqn,info] = fourCurl3P2(node,elem,bdFlag,pde,option);
    % compute error
    uI = edgeinterpolate2(pde.exactu,node,eqn.edge,eqn.face,eqn.face2edge);
    wI = edgeinterpolate2(pde.curlcurlu,node,eqn.edge,eqn.face,eqn.face2edge);
    errwL2(k) = getL2error3NE2(node,elem,pde.curlcurlu,w);
    erruL2(k) = getL2error3NE2(node,elem,pde.exactu,u);   
    errwwhHcurl(k) = getHcurlerror3NE2(node,elem,pde.curlcurlcurlu,w);
    erruuhHcurl(k) = getHcurlerror3NE2(node,elem,pde.curlu,u); 
    errwIwhL2(k) = sqrt((w-wI)'*eqn.M*(w-wI));
    erruIuhL2(k) = sqrt((u-uI)'*eqn.M*(u-uI)); 
    errwIwhHcurl(k) = sqrt((w-wI)'*eqn.A*(w-wI));
    erruIuhHcurl(k) = sqrt((u-uI)'*eqn.A*(u-uI));      
    % record information
    N(k) = length(u)+length(w);
    h(k) = 1./(size(node,1)^(1/3)-1);            
    assembleTime(k) = info.assembleTime;
    solverTime(k) = info.solverTime;
end

%% Plot convergence rates
subplot(2,2,1)
showrateh2(h, errwIwhL2(1:maxIt), 1 , 'k-+', '||w_I-w_h||',...
           h, erruIuhL2(1:maxIt), 1 , 'r-*', '||u_I-u_h||' );
subplot(2,2,2)
showrateh2(h, errwL2(1:maxIt), 1 , 'k-+', '||w-w_h||',...
           h, erruL2(1:maxIt), 1 , 'r-*', '||u-u_h||' );                                     
subplot(2,2,3)               
showrateh2(h, errwIwhHcurl(1:maxIt), 1 , 'b-+', '||w_I-w_h||_1',...
           h, erruIuhHcurl(1:maxIt), 1 , 'g-*', '||u_I-u_h||_1' );
subplot(2,2,4)               
showrateh2(h, errwwhHcurl(1:maxIt), 1 , 'b-+', '||w-w_h||_1',...
           h, erruuhHcurl(1:maxIt), 1 , 'g-*', '||u-u_h||_1' );
                                                                           
%% Output
err = struct('N', N, 'h', h, ...
             'wL2',errwL2(1:maxIt),'uL2',erruL2(1:maxIt),...
             'wIwhL2',errwIwhL2(1:maxIt), 'uIuhL2',erruIuhL2(1:maxIt),...
             'wIwhHcurl',errwIwhHcurl(1:maxIt),'uIuhHcurl',erruIuhHcurl(1:maxIt),...
             'wwhHcurl',errwwhHcurl(1:maxIt),'uuhHcurl',erruuhHcurl(1:maxIt));
time = struct('N',N,'assemble',assembleTime(1:maxIt),'solver',solverTime(1:maxIt));
       
%% Display error on screen
disp('Table: Error')
colname = {'#Dof','h','||w_I-w_h||','||w-w_h||','||w_I-w_h||_1','||w-w_h||_1'};
disptable(colname,err.N,[],err.h,'%0.3e',err.wIwhL2,'%0.5e',err.wL2,'%0.5e',...
              err.wIwhHcurl,'%0.5e',err.wwhHcurl,'%0.5e');

colname = {'#Dof','h','||u_I-u_h||','||u-u_h||','||u_I-u_h||_1','||u-u_h||_1'};
disptable(colname,err.N,[],err.h,'%0.3e',err.uIuhL2,'%0.5e',err.uL2,'%0.5e',...
              err.uIuhHcurl,'%0.5e',err.uuhHcurl,'%0.5e');    
          
%     display('Table: CPU time')
colname = {'#Dof','Assemble','Solve'};
disptable(colname,time.N,[],time.assemble,'%0.2e',time.solver,'%0.2e');
          