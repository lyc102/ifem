function t2v = cell2node(elem)
%% CELL2NODE computes the sparse matrix representing elem for polygonal mesh
% elem{i}: i-th elem vertex indices, elem is a cell array
% Reference: A Simple Virtual Element-based Flux Recovery on Quadtree, 
%            https://arxiv.org/abs/2006.05585, ERA 2021
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

elemVertNum = cellfun('length',elem);
minNv = min(elemVertNum);
maxNv = max(elemVertNum);
NT = size(elem,1);
N = max([elem{:}]); % same with N = max(cellfun(@max, elem));
t2v = sparse(NT,N);
for nV = minNv:maxNv
    isNv = (elemVertNum == nV);
    if ~(any(isNv)); continue; end
    ixNv = find(isNv);
    elemNv = cell2mat(elem(isNv));
    t2v = t2v + sparse(repmat(ixNv, [nV, 1]), elemNv(:), 1, NT,N);
end