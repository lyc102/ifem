function [Ac,Pro,Res] = coarsenAMGrs(A,theta)

% load lakemesh % For debug

if nargin<2, theta = 0.025; end

N = size(A,1);
%% Generate strong connectness matrix
maxaij = min(A); % max of negative of off-diagonals
% if min(maxaij) > 0 % all off-diagonals are positive
%     isF = true(N,1); isC = [];
% end
D = spdiags(1./abs(maxaij'),0,N,N);
Am = D*A;

%% Delete weak connectness
[im,jm,sm] = find(Am); 
idx = (-sm > theta); % delete weak connectness and diagonal
As = sparse(im(idx),jm(idx),1,N,N); % asymmetric matrix for directed graph
Am = sparse(im(idx),jm(idx),sm(idx),N,N);
Ass = (As + As')/2;  % symmetric matrix for the undirected graph

% %% Normalize A
% D = spdiags(1./sqrt(diag(A)),0,N,N);
% Am = D*A*D;
% 
% %% Delete weak connectness
% [im,jm,sm] = find(Am); 
% idx = (-sm > theta); % delete weak connectness and diagonal
% As = sparse(im(idx),jm(idx),1,N,N);
% % As = sparse(im(idx),jm(idx),sm(idx));
% Ass = As;

%% Put isolated nodes to F set
isF = false(N,1);     % F: fine node
degIn = sum(As);      % number of vertices strong connected to
degIn = full(degIn');
isF(degIn == 0) = true; % isolated nodes are fine nodes

%% Find an approximate maximal independent set and put to C set
isC = false(N,1);       % C: coarse node
U = (1:N)';             % U: undecided node
degFin = zeros(N,1);    % number of F nodes strong connected to
% debug
% close all;
% scrsz = get(0,'ScreenSize');
% figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
% showmesh(node,elem); 
while sum(isC) < N/2 && length(U)>20
    % Randomly pick up around half vertex of undecided nodes
    % the bigger degree is, the higher possibility
    isS = false(N,1);  % S: selected set, changing in the coarsening
    degInAll = degIn + degFin;
    isS((rand(N,1) < 0.85*degInAll/mean(degInAll)) & (degInAll>0)) = true; 
    S = find(isS); % transfer index for S to U
%     S = U(SinUidx);
    % debug
    % findnode(node,S,'Color','g');

    % Delete one vertex from connected pair in S
    [i,j] = find(triu(Ass(S,S),1));    % i,j: index for S
%     isS(S(j)) = false;
    idx = degInAll(S(i)) >= degInAll(S(j)); % compare degree of vertices
    isS(S(j(idx))) = false; % remove vertices with smaller degree 
    isS(S(i(~idx))) = false;
    isC(isS) = true;
    % findnode(node,isS,'Color','r');

%     fprintf('\n Number of independent points: %6.0u\n',sum(isS));

    % Add neighboring nodes of C points into F set
    [i,j] = find(Ass(:,isC)); %#ok<*NASGU>
    isF(i) = true;
    U = find(~(isF | isC)); % the rest is undecided
    
    % update degree (weight) of each node
    degIn(isF | isC) = 0;    % 
    degFin = zeros(N,1);
    degFin(U) = full(sum(As(isF,U))); % strong connectness with F points
    
    if length(U) <= 20 % add left nodes into C nodes
        isC(U) = true;
        U = [];       % to exit the while loop;
    end
    % debug
%     showmesh(node,elem); 
    % findnode(node,isF,'Color','y')
    % findnode(node,U,'Color','m')
    % findnode(node,isC,'Color','r');    
end
fprintf('Number of coarse nodes: %6.0u\n',sum(isC));

%% Coarse grid and fine grid index map
allNode = (1:N)';
fineNode = allNode(~isC);
Nf = length(fineNode);
Nc = N-length(fineNode);
coarseNode = (1:Nc)';
coarse2fine = find(isC);
fine2coarse = zeros(N,1);
fine2coarse(isC) = coarseNode;
ip = coarse2fine; 
jp = coarseNode;
sp = ones(Nc,1);
% save indexmap coarse2fine fine2coarse

%% Construction prolongation operator
Afc = Am(fineNode,coarse2fine);
% Afc = transpose(Am(coarse2fine,fineNode));
Dsum = sum(Afc,2);
% Dsum = sparse(1:Nf,1:Nf,1./Dsum);
Dsum = spdiags(1./Dsum, 0, Nf, Nf);
[ti,tj,tw] = find(Dsum*Afc);
ip = [ip; fineNode(ti)];  % fine node index
jp = [jp; tj]; % coarse node index
sp = [sp; tw];
Pro = sparse(ip,jp,sp,N,Nc);
Res = Pro';

%% Form coarse grid matrix
Ac = Res*A*Pro;