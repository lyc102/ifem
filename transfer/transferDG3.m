function [d2c, c2d] = transferDG3(elem2dof)
% TRANSFERDG3 constructs the sparse transfer matrix from a discrete data
% structure (DG, double-valued face-hybrid variable) to a continuous data
% structure (continuous nodal DoFs, singled-valued face-hybrid). The discrete
% DoF is assumed to be 
%
%     elem2discreteDof = reshape(1:numDofK*NT, numDofK, NT)'
% 
% [d2c, c2d] = transferDG3(elem) computes the transfer matrix from the
% discrete nodal dofs to continuous or vice versa, suppose one assembles a
% bilinear form matrix B_{ij}:= B(u_j,v_i) using a discrete space, then the 
% continuous one can be simply done by 
%
%       B = c2d'*B*c2d
% 
% [d2c, c2d] = transferDG3(elem2face) compute the transfer matrix
% from the discrete face dofs to continuous.
%
% Remark: note a more readable d2c for face code is as follows but less concise
%       face2elemIdx1 = (face2elem(:,1)-1)*4+face2elem(:,3);
%       face2elemIdx2 = (face2elem(:,2)-1)*4+face2elem(:,4);
%       isIntFace = (face2elem(:,1) ~= face2elem(:,2));
%       d2c = sparse((1:NF)', face2elemIdx1, ones(NF, 1), NF, 4*NT);
%       d2c = d2c + sparse((1:NF)', face2elemIdx2, isIntFace, NF, 4*NT);
%
% See also Poisson3Hybrid, coarsen
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% 
NT = size(elem2dof,1);
Ndof = max(elem2dof(:));
numDofK = size(elem2dof,2);
elem2disdof = reshape(1:numDofK*NT, numDofK, NT)';
c2d = sparse(elem2disdof(:), elem2dof(:), 1, numDofK*NT, Ndof);
d2c = c2d';
end