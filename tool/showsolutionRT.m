function showsolutionRT(node,elem,u,varargin)
%% SHOWSOLUTIONRT plots a RT function u on a triangular mesh in 2-D.
%
%  [node, elem] = squaremesh([0,1,0,1],1/2^6);
%  pde = sincosdata;
%  sigmaI = faceinterpolate(pde.Du,node,elem,'RT0');
%  showsolutionRT(node,elem,sigmaI);
%
%  See also showsolution
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

NT = size(elem,1);

%% Generate a big triangulation
elemnew = reshape(1:3*NT,NT,3);
nodenew = node(elem(:),:);

%% Evaluate piecewise linear function
elem = sortelem(elem);
elem2edge = dofedge(elem);
Clambda = curlbasis(node,elem);
unew(1:NT,:) = [u(elem2edge(:,3)) u(elem2edge(:,3))].*Clambda(:,:,2) ...
             + [u(elem2edge(:,2)) u(elem2edge(:,2))].*Clambda(:,:,3);
unew(NT+(1:NT),:) = - [u(elem2edge(:,3)) u(elem2edge(:,3))].*Clambda(:,:,1) ...
                    + [u(elem2edge(:,1)) u(elem2edge(:,1))].*Clambda(:,:,3);
unew(2*NT+(1:NT),:) = - [u(elem2edge(:,2)) u(elem2edge(:,2))].*Clambda(:,:,1) ...
                      - [u(elem2edge(:,1)) u(elem2edge(:,1))].*Clambda(:,:,2);

%% Plot the solution
subplot(1,2,1); showsolution(nodenew,elemnew,unew(:,1),varargin{:});
subplot(1,2,2); showsolution(nodenew,elemnew,unew(:,2),varargin{:});