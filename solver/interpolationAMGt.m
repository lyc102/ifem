function [Pro,Res] = interpolationAMGt(A,isC)   
%% INTERPOLATIONAMGT construct prolongation and restriction matrices
%
% [Pro,Res] = INTERPOLATIONAMGT(A,isC) construct prolongation and
% restriction matrices use matrix-dependent interpolation. Each fine nodes
% use at most two coarse nodes.
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
Afc = A(fineNode,coarse2fine);     % matrix-dependent interpolation
[i,j,s] = find(Afc);
i2j1(i(1:end)) = j(1:end);
i2s1(i(1:end)) = s(1:end);
i2j2(i(end:-1:1)) = j(end:-1:1);
i2s2(i(end:-1:1)) = s(end:-1:1);
Dsum = i2s1' + i2s2';
idx = (Dsum ~= 0);
Dsum = Dsum(idx);
i2j1 = i2j1(idx)';  i2s1 = i2s1(idx)';
i2j2 = i2j2(idx)';  i2s2 = i2s2(idx)';
% Dsum = sum(Afc,2);
% Dsum = spdiags(1./Dsum, 0, Nf, Nf); % normalize to preserve constant
% [ti,tj,tw] = find(Dsum*Afc);
% ip = [ip; fineNode(ti)];  % fine node index
% jp = [jp; tj];            % coarse node index
% sp = [sp; tw];            % weight
ip = [ip; fineNode(idx); fineNode(idx)];  % fine node index
jp = [jp; i2j1; i2j2];            % coarse node index
sp = [sp; i2s1./Dsum; i2s2./Dsum];            % weight
Pro = sparse(ip,jp,sp,N,Nc);
Res = Pro';