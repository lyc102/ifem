function [L,p,Ac] = acholold(A)

N = size(A,1);
%% Multilvel factorization
% level = ceil(log(N)/log(2))+1;
level = 3;
R = cell(level,1);
isM = cell(level,1);
Nl = zeros(level,1);
Nf = zeros(level,1);
for k = 1:level-1
    Nl(k) = size(A,1);    
    [isM{k},As] = mis(A);
    Nf(k) = sum(isM{k});
    [R{k},A] = afold(A,As,isM{k});
end
% coarest level
Nl(level) = size(A,1);
Nf(level) = Nl(level);
isM{level} = true(Nl(level),1);
R{level} = speye(Nl(level),Nl(level));
% R{level} = spdiags(sqrt(diag(A)),0,Nl(level),Nl(level));
Ac = A; % coarest grid matrix
% R{level} = chol(A);
% opts.type = 'nofill'; opts.michol = 'on';
% R{level} = ichol(A,opts);

%% Ordering/Permutation for the lower part
pl = cell(level,1);
pl{level-1} = 1:Nl(level);
for k = level-2:-1:1
    Ck = (1:Nl(k+1))';     % index for coarse nodes in level k
    Ck1 = Ck(isM{k+1});    % fine nodes in level k+1  goes first
    Ck2 = Ck(~isM{k+1});
    pl{k} = [Ck1; Ck2(pl{k+1})];
           % Ck(isM{k+1})  
           % Ck(P{k+1})    map ordering of level k+1 to level k
    R{k}(:,Ck+Nf(k)) = R{k}(:,pl{k}+Nf(k));
end
coarseNode = find(~isM{1});
p = [find(isM{1}); coarseNode(pl{1})];

%% Merge multilevel R to one lower triangular matrix
L = sparse(N,N);
L(1:N,1:Nf(1)) = transpose(R{1});
col = cumsum(Nf);
for k = 2:level
    L(col(k-1)+1:N,col(k-1)+1:col(k)) = transpose(R{k});
end

%% Merge multilevel factorization into a big one
% record permuation of nodes
% p = zeros(N,1);
cpointer = 1;
pc = (1:N)'; % global index;
for k = 1:level-1    
    endpointer = cpointer + sum(isM{k})-1;
    p(cpointer:endpointer) = pc(isM{k});
    pc = pc(~isM{k});     % global indx of coarse nodes
    cpointer = endpointer+1;
end
p(cpointer:N) = pc;

% generate a tempor i,j,s for debug
i = cell(level,1);
j = cell(level,1);
s = cell(level,1);
for k = 1:level
    tempL = transpose(R{k});
    [i{k},j{k},s{k}] = find(tempL);
end

% shift and merge index in from coarse level to fine level
for k = level:-1:2
    Nf = sum(isM{k-1});
    % shift the index in coarse level
    i{k} = i{k} + Nf;
    j{k} = j{k} + Nf;
    % add to fine level, i.e. merge two levels
    i{k-1} = [i{k-1}; i{k}];
    j{k-1} = [j{k-1}; j{k}];
    s{k-1} = [s{k-1}; s{k}];
end
% merge two level L into one
L = sparse(i{1},j{1},s{1},N,N);