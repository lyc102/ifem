function [eta,Du] = estimaterecovery3(node,elem,u)
% ESTIMATERECOVERY3 recovery type error estimator in 3-D.
%  
% eta = estimaterecovery3(node,elem,u) computes an error estimator eta by
% recovey second derivative of the gradient of a finite element function u. 
%
% [eta, Du] = estimaterecovery3(node,elem,u) also returns the recovered
% derivative Du which is in P1 finite element space.
%
% By interpolation error estimate $|u-uI|_{1,2,T}\leq Ch_T|u|_{2,2,T}$.
% Therefore we recovery an approximation of second derivatives of u and
% compute elementwise H2 norm. We use the weighted averaging recovery
% scheme with volume weight.
%
% See also recovery3, estimaterecovery, cube
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

[Dlambda,volume] = gradbasis3(node,elem);
Du = gradu3(node,elem,u,Dlambda);
Du = recovery3(node,elem,Du,volume);
DDu(:,1:3) = gradu3(node,elem,Du(:,1),Dlambda);
DDu(:,4:6) = gradu3(node,elem,Du(:,2),Dlambda);
DDu(:,7:9) = gradu3(node,elem,Du(:,2),Dlambda);
h = volume.^(1/3);
eta = h.*sqrt(volume.*sum(DDu.^2,2));
%% TODO a test example
% TODO add diffusion coefficient