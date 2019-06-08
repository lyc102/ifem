function [node2agg,As] = coarsenAMGa(A,theta)

%% Parameters
if ~exist('theta','var'), theta = 0.025; end
N = size(A,1);
isC = false(N,1);       % C: coarse node
N0 = min(sqrt(N),25);   % number of the coarest nodes

%% Initialize output
node2agg = zeros(N,1);
agg2node = zeros(N,1);

%% Generate strong connectness matrix
Dinv = spdiags(1./sqrt(diag(A)),0,N,N);
Am = Dinv*A*Dinv;       % normalize diagonal of A
[im,jm,sm] = find(Am); 
idx = (-sm > theta);    % delete weakly connect off-diagonal and diagonal
As = sparse(im(idx),jm(idx),sm(idx),N,N); % matrix for strong connectness
As = As + speye(N);    % add diagonal
As1 = spones(As);      % graph of As
As2 = triu(As1*As1,1); % edges of the graph corresponding to As^2

%% Compute degree of vertex
deg = full(transpose(sum(As1))); % number of strongly connected neighbors
% the case for mass matrix requires different treatment: To do
if sum(deg>0) < 0.25*sqrt(N)  % too few connected nodes e.g. A is mass matrix
    isC(ceil(rand(N0,1)*N)) = true; % randomly chose N0 nodes
    agg2node = find(isC);
    node2agg(isC) = (1:length(agg2node))';
    return                    % smoother is a good preconditioner
end           
idx = (deg>0);
deg(idx) = deg(idx) + 0.1*rand(sum(idx),1); % break the equal degree case

%% Find an approximate maximal independent set and put to C set
isC = false(N,1);       % C: coarse node
isF = false(N,1);       % F: fine node
isU = true(N,1);        % U: undecide node
isS = true(N,1);        % S: selected node
isF(deg == 0) = true;   % isolated nodes are added into F set
aggN = 0;               % aggregrate number
while aggN < N/2 && sum(isS) >N0 
    % Mark all undecided nodes
    isS = false(N,1);  % S: selected set, changing in the coarsening
    isS(deg>0) = true;
    S = find(isS); 
    
    % Find marked nodes with local maximum degree
    [i,j] = find(As2(S,S));   % i,j and i<j: edges of subgraph S
    idx = deg(S(i)) >= deg(S(j));     % compare degree of vertices
    isS(S(j(idx))) = false;  % remove vertices with smaller degree 
    isS(S(i(~idx))) = false; % if degrees are equal, keep the nodes with smaller index 
    isC(isS) = true;                % set selected nodes as coarse nodes
    
    % Add new agg
    newC = find(isS);
    newAgg = aggN+(1:length(newC));
    aggN = aggN + length(newC);
    node2agg(newC) = newAgg;  
    agg2node(newAgg) = newC;
    
    % Remove coarse nodes and add neighboring nodes to the aggregate
    U = find(isU);
    [i,j] = find(As(isU,newC)); %#ok<*NASGU> use original connectivity
    isF(U(i)) = true;        % neighbor of C nodes are F nodes
    isU = ~(isF | isC);      % U: undecided set
    node2agg(U(i)) = node2agg(newC(j)); % add neighbors into the same aggregate
    deg(newC) = 0;           % remove new selected coarse grid nodes
    deg(U(i)) = 0;           % remove neighbors of new selected coarse grid nodes
    U = find(isU);
    [i,j] = find(As(U,isF)); % find neighbor of fine nodes
    deg(U(i)) = 0;              % remove neighbors of existing agg
    
%     showmesh(node,elem); 
%     showagg(node,node2agg,agg2node,As);
end
agg2node = agg2node(1:max(node2agg));

%% Add left vertices to existing agg
while any(isU)
    U = find(isU);           % undecided nodes
    [i,j] = find(As(:,isU)); %#ok<*NASGU> neighboring nodes of U
    % j: undecided vertices; i: neighbor of j
    neighborAgg = node2agg(i);% agg number of neighbors
    idx = (neighborAgg > 0);  % a interior nodes could be left
    [nAgg,neighborAgg] = max(sparse(neighborAgg(idx),j(idx),1));   
    % a undecided node could be next to several aggregrates. find the
    % aggregate to which this node containing maximal edges.
    isbdU = (nAgg > 0);      % find undecided nodes next to 
    bdU = U(isbdU);          % the boundary of aggregates
    node2agg(bdU) = neighborAgg(isbdU); % assgin the aggregate number 
    % remove bdU from U and add to F
    isF(bdU) = true;
    isU(bdU) = false;
    
%     findnode(node,bdU,'noindex','Color','m');
%     showagg(node,node2agg,agg2node,As);
end