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
HB(:,1) = fineNode;
% Nf = length(fineNode);  
Nc = N-length(fineNode);
coarseNode = (1:Nc)';          % coarse node index
coarseNodeFineIdx = find(isC); % coarse node index in the fine grid
% fine2coarse(isC) = coarseNode;

%% Construct matrix-dependent prolongation and restriction operator
Afc = A(fineNode,coarseNodeFineIdx);     % sub-graph between F and C
[i,j,s] = find(Afc);                     % edges between [F C]
HB(i(end:-1:1),3) = j(end:-1:1);         % random chose two nodes in C
HB(i(1:end),2) = j(1:end);               % as parents of a node in F
w(i(end:-1:1),2) = s(end:-1:1);          
w(i(1:end),1) = s(1:end);                % use corresponding aij as weight 
Dsum = sum(w,2);     
idx = (Dsum ~= 0);                       % consider non isolated vertices
w = w(idx,:)./[Dsum(idx), Dsum(idx)];    % normalize the weight
HB = HB(idx,:);
ip = [coarseNodeFineIdx; HB(:,1); HB(:,1)];  % fine nodes index
jp = [coarseNode;        HB(:,2); HB(:,3)];  % coarse nodes index
sp = [ones(Nc,1);         w(:,1); w(:,2)];   % weight
Pro = sparse(ip,jp,sp,N,Nc);
Res = Pro';