function [c2d, d2c] = transferDG3(elem2dof)
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
% [c2d, d2c] = transferDG3(elem2face) compute the transfer matrix
% from the discrete face dofs to continuous.
%
% Remark: note a more readable d2c for face code is as follows but less concise
%       face2elemIdx1 = (face2elem(:,1)-1)*4+face2elem(:,3);
%       face2elemIdx2 = (face2elem(:,2)-1)*4+face2elem(:,4);
%       isIntFace = (face2elem(:,1) ~= face2elem(:,2));
%       d2c = sparse((1:NF)', face2elemIdx1, ones(NF, 1), NF, 4*NT);
%       d2c = d2c + sparse((1:NF)', face2elemIdx2, isIntFace, NF, 4*NT);
%
% Example 1: construct the transfer matrix from discrete face quadratic Nedelec dof 
% to continuous face quadratic Nedelec dof, given elem2face(t, i) represents the
% global indexing of the t-th element's i-th face
% 
%       elem2dof = dof3expand(elem2face, 12); % 12 dofs per face for Nd2
%       c2d = transferDG3(elem2dof);
% 
% Example 2: compute the singled valued face dof from the double valued face dof
% uh is a variable of size (4*NT,1) that defines the double valued face dof, then
% uhat is {{uh}}_F is the single valued face dof.
% 
%       c2d = transferDG3(elem2face);
%       uhat = 0.5*c2d'*uh;
%
% See also Poisson3Hybrid, dof3expand, coarsen
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