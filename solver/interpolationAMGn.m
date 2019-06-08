function [Pro,Res] = interpolationAMGn(A,isC)   
%% INTERPOLATIONAMGS construct prolongation and restriction matrices
%
% [Pro,Res] = interpolatoinAMGs(A,isC) construct prolongation and
% restriction matrices use standard matrix-dependent interpolation. 

% In the input, A is a SPD matrix and isC is a logical array to indicate
% nodes in coarse matrix. In the output Pro and Res are prolongation and
% restriction matrices satisfying Res = Pro'.
%
% Example
%   load lakemesh
%   A = assemblematrix(node,elem);
%   [isC,As] = coarsenAMGc(A);
%   [Ac,Pro,Res] = interpolatoinAMGs(As,isC);
%
% See also: coarsenAMGc, amg
% 
% New interpolation added by Xiaozhe Hu. 12/10/2011.


N = size(A,1);

%% Index map between coarse grid and fine grid
allNode = (1:N)';          
fineNode = allNode(~isC);
Nf = length(fineNode);  
Nc = N-length(fineNode);
coarseNode = (1:Nc)';    % coarse node index
coarse2fine = find(isC); % coarse node index in the fine grid
% fine2coarse(isC) = coarseNode;
ip = coarse2fine;  % coarse node index in the fine grid
jp = coarseNode;   % coarse node index
sp = ones(Nc,1);   % identity matrix for coarse nodes in the fine grid

%% Construct prolongation and restriction operator
Afc = A(fineNode, coarse2fine);
Dsum = spdiags(A(fineNode,fineNode),0);
alpha = (sum(A(fineNode,:),2) - Dsum)./sum(Afc,2);
Dsum = spdiags(((-1).*alpha)./Dsum, 0, Nf, Nf); 
[ti,tj,tw] = find(Dsum*Afc);
ip = [ip; fineNode(ti)];  % fine node index
jp = [jp; tj];            % coarse node index
sp = [sp; tw];            % weight
Pro = sparse(ip,jp,sp,N,Nc);
Res = Pro';