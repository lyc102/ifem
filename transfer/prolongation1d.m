function Ix = prolongation1d(level)
%% PROLONGATION1D 
% construct a sparse matrix representation of the prolongation operator of
% one dimensional vector. it maps a vector 2^level to 2^(level+1). The
% default is zero boundary condition.

N = 2^level+1; % number of points
Nf = 2^(level+1) + 1; % number of points in the fine grid
j = 2:N-1;     % interiori points
jf = 2*j-1; 
ii = [jf jf-1 jf+1];
jj = [j j j];
ss = [ones(1,N-2) 0.5*ones(1,N-2) 0.5*ones(1,N-2)];
Ix = sparse(ii,jj,ss,Nf,N);