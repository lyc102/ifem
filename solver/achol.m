function [L,p,Ac] = achol(A)
%% ACHOL approximated Cholosiky factorization
%
% L: lower triangular matrix
% p: permutation s.t. L is an approximated factorization of A(p,p)
% p(1:Nf) is the index of fine nodes

N = size(A,1);
p = zeros(N,1);

%% Multilvel factorization
level = ceil(log(N)/log(2))+1;
% level = 3;
i = cell(level,1);
j = cell(level,1);
s = cell(level,1);
isM = cell(level,1);
cpointer = 1;
pc = (1:N)'; % global index;
for k = 1:level-1
    % find Af and Afc
    [isM{k},As] = mis(A);
    [i{k},j{k},s{k},A] = af(A,As,isM{k});
    % map the index to global index of nodes
    i{k} = pc(i{k});
    j{k} = pc(j{k});
    % record permutation
    endpointer = cpointer + sum(isM{k})-1;
    p(cpointer:endpointer) = pc(isM{k});
    pc = pc(~isM{k});     % global indx of coarse nodes
    cpointer = endpointer+1;  
    % exit the multilevel decomposition if A is empty
    if isempty(A)
        level = k;
        break
    end
end
p(cpointer:N) = pc;
% coarest level
Ac = A; % coarest grid matrix
if ~isempty(Ac)
    Nc = size(A,1);
    i{level} = pc((1:Nc)');
    j{level} = pc((1:Nc)');
    s{level} = ones(Nc,1);
end
% Lc = ichol(Ac);
% [ic,jc,sc] = find(Lc);
% i{level} = pc(ic);
% j{level} = pc(jc);
% s{level} = sc;

%% Merge multilevel factorization into a big one
% map the index to permuated one 
pinv(p) = (1:N)';
ii = pinv(cell2mat(i));
jj = pinv(cell2mat(j));
ss = cell2mat(s);
% generate L
L = sparse(ii,jj,ss,N,N);
% R = L';