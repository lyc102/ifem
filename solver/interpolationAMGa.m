function [Pro,Res] = interpolationAMGa(A,isC)   
%% INTERPOLATIONAMGA construct prolongation and restriction matrices
%
% [Pro,Res] = INTERPOLATIONAMGA(A,isC) construct prolongation and
% restriction matrices use matrix-dependent interpolation. Each fine nodes
% use only one coarse node.
%
% In the input, A is a SPD matrix and isC is a logical array to indicate
% nodes in coarse matrix. In the output Pro and Res are prolongation and
% restriction matrices satisfying Res = Pro'.
%
% The submatrix A_{cf} is used to construct the interpolation of values on
% fine nodes from that of coarse nodes. The weight is normalized to
% preserve the constant.
%
% Example
%   load lakemesh
%   A = assemblematrix(node,elem);
%   [isC,As] = coarsenAMGc(A);
%   [Pro,Res] = interpolation(As,isC);
%
% See also: coarsenAMGc, amg
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

N = size(A,1);
%% Index map between coarse grid and fine grid
allNode = (1:N)';          
fineNode = allNode(~isC);
% Nf = length(fineNode);  
Nc = N-length(fineNode);
coarseNode = (1:Nc)';    % coarse node index
coarse2fine = find(isC); % coarse node index in the fine grid
% fine2coarse(isC) = coarseNode;
ip = coarse2fine;  % coarse node index in the fine grid
jp = coarseNode;   % coarse node index
sp = ones(Nc,1);   % identity matrix for coarse nodes in the fine grid

%% Construct prolongation and restriction operator
Afc = -A(fineNode,coarse2fine);     % matrix-dependent interpolation
[Dsum,j] = max(Afc,[],2);
idx = (Dsum ~= 0);
ip = [ip; fineNode(idx)];  % fine node index
jp = [jp; j(idx)];         % coarse node index
sp = [sp; ones(length(j(idx)),1)];    % weight
Pro = sparse(ip,jp,sp,N,Nc);
Res = Pro';