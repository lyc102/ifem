function [Pro,Res] = interpolationAMGsa(A,node2agg,omega,smoothingstep)   
%% INTERPOLATIONAMGSA construct prolongation using smoothed aggregation
%
% [Pro,Res] = INTERPOLATIONAMGSA(A,node2agg) construct prolongation and
% restriction matrices using smoothed aggregation.
%
% In the input, A is a SPD matrix and node2agg assigns an aggregate for
% each node in A. In the output Pro and Res are prolongation and
% restriction matrices satisfying Res = Pro'.
%
% The prolongation operator is first set to be piecewise constant and then
% using matrix A to smooth out this simple prolongation using weighted
% Jacobi iteration for several steps. The default choice of the weight is
% 0.35 and the smoothingstep is 2. Larger smoothing steps will result a
% denser prolongation operator.
%
% Example
%   load lakemesh
%   A = assemblematrix(node,elem);
%   [node2agg,As] = coarsenAMGa(A);
%   [Pro,Res] = interpolationAMGsa(As,node2agg);
%
% See also: coarsenAMGc, amg
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if ~exist('omega','var'), omega = 0.35; end
if ~exist('smoothingstep','var'), smoothingstep = 2; end

%% A simple prolongation
N = size(A,1);
Nc = max(node2agg);
idx = find(node2agg~=0);
Pro = sparse(idx,node2agg(idx),1,N,Nc);

%% Smooth the piecewise constant prolongation
for k = 1:smoothingstep
    Pro = Pro - omega*(A*Pro);
end

%% Normalize the prolongation such that the row sum is one
rowsum = sum(Pro,2);
D = spdiags(1./rowsum,0,N,N);
Pro = D*Pro;
Res = Pro';
