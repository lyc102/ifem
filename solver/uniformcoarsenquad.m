function [elem,HB] = uniformcoarsenquad(elem)
%% UNIFORMCOARSENQUAD uniform coarsening of uniform quad mesh
%
% [elem,HB] = UNIFORMCOARSENQUAD(elem) remove grid points added by uniform
% refinement. See the illustration below:
%
%
%   See also: uniformcoarsen, coarsen, bisect, uniformcoarsen3, mg
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

NT = size(elem,1);
N = max(elem(:));
d = size(elem,2);
if d == 3  % triangulation
    NT = NT/2;
end
HB = [];

%% Find ni and nj
b = N - NT + 1;
c = N;
delta = sqrt(b^2-4*c); 
ni = (b+delta)/2;
nj = (b-delta)/2;
if ~isequal(mod(ni+1,2),0) || ~isequal(mod(nj+1,2),0)
    display('Not from refinement of uniform quad mesh');
    return
end
if d == 3 % triangulation
    if ~isequal(elem(1,1),ni+1)
        nj = ni;
        ni = elem(1,1) - elem(1,2);
    end    
elseif d == 4
    if ~isequal(elem(1,2),ni+1)
        nj = ni;
        ni = elem(1,2) - elem(1,1);
    end        
end

%% Coarse grids
nic = (ni+1)/2;
njc = (nj+1)/2;
Nc = nic*njc;
nodecidx = reshape(1:Nc,nic,njc);
t2nidxMap = nodecidx(1:nic-1,1:njc-1);
k = t2nidxMap(:);
elem = [k k+nic k+nic+1 k+1];

%% Record HB
HB = zeros(N,3);
nodeidx = reshape(1:N,ni,nj);
% fine nodes on vertical lines
i = 2:2:ni-1;
j = 1:2:nj;
k = nodeidx(i,j);
HB(k,1) = k(:);
HB(k,2) = k(:) - 1;
HB(k,3) = k(:) + 1;
% fine nodes on horizontal lines
i = 1:2:ni;
j = 2:2:nj-1;
k = nodeidx(i,j);
HB(k,1) = k(:);
HB(k,2) = k(:) - ni;
HB(k,3) = k(:) + ni;
% fine nodes in the center
i = 2:2:ni-1;
j = 2:2:nj-1;
k = nodeidx(i,j);
HB(k,1) = k(:);
HB(k,2) = k(:) - ni-1;
HB(k,3) = k(:) + ni+1;
HB(HB(:,2) == 0,:) = [];
% shift the index into coarse grid
i = 1:2:ni;
j = 1:2:nj;
k = nodeidx(i,j);
indexMap = zeros(N,1);
indexMap(k(:)) = 1:Nc;
HB(:,2:3) = indexMap(HB(:,2:3));