function Aij = spsub2ind(A,i,j)
%% SPSUB2IND find values of a sparse matrix 
%
% Aij = SPSUB2IND(A,i,j) extract values of the sparse matrix A with
% subscrpts (i,j) where (i,j) could be arrays with length > 1. The intutive
% way A(i,j) will give a sub-matrix.
% 
% Example
%
%    n = 5;
%    e = ones(n,1);
%    A = spdiags([e -2*e e], -1:1, n, n);
%    A(3,2) = 5;
%    i = [1 2 3];
%    j = [1 3 2];
%    display(full(A));
%    display(spsub2ind(A,i,j));
%
% See also sub2ind, ind2sub
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

[M,N] = size(A);
if size(i,1) == 1, i = i'; end    % change row vector to column vector
if size(j,1) == 1, j = j'; end
idx = [i,j];
[sortedidx,I] = sortrows(idx,2);  %#ok<*ASGLU> % sort the subscript by j
sppattern = sparse(i,j,1,M,N);    % sparse pattern using i,j subscripts
[ii,jj,zz] = find(sppattern.*A);  % sppattern.*A will keep only (i,j) values of A
% But for some (i1,j1), A(i1,j1) = 0 and thus not in the output of find.
% We generate a sparse matrix using (i,j) and (ii,jj) subscripts. In newA,
% newA(i1,j1) = max(abs(zz)) + 1 is nonzero and still keep other nonzeros.
newA = sparse([i; ii],[j; jj],[max(abs(zz))+ ones(length(i),1); zz]);
[si,sj,sortedAij] = find(newA);   % [si,sj] is sorted by sj
Aij(I) = sortedAij - 1 - max(abs(zz));