function [eta,Du] = estimaterecovery(node,elem,u)
%% ESTIMATERECOVERY recovery type error estimator.
%  
% eta = estimaterecovery(node,elem,u) computes an error estimator eta by
% recovery second derivative of the gradient of a finite element function u. 
%
% [eta, Du] = estimaterecovery(node,elem,u) also returns the recovered
% derivative Du which is in P1 finite element space.
%
% By interpolation error estimate $|u-uI|_{1,2}\leq C|u|_{2,1}$. Therefore
% we recovery an approximation of second derivatives of u and compute L1
% norm. We use the weighted averaging recovery scheme with area weight.
%
% See also recovery, estimaterecovery3, Lshape, crack
% 
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

[Dlambda,area] = gradbasis(node,elem);
Du = gradu(node,elem,u,Dlambda);
Du = recovery(node,elem,Du,area);
DDu(:,1:2) = gradu(node,elem,Du(:,1),Dlambda);
DDu(:,3:4) = gradu(node,elem,Du(:,2),Dlambda);
eta = area.*sum(abs(DDu),2);
%% TODO add diffusion coefficient