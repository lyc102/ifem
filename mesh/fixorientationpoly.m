function elem = fixorientationpoly(node,elem,varargin)
%% FIXORIENTATIONPOLY set all polygon clockwise or counter-clockwise
% 
%   elem = fixorientationpoly(node,elem) reorder the elem to orient
%   clockwise similar to the quad element set up. 
%   elem{:} has to be a cell array. Default is counter-clockwise order.
%   Reference: A simple virtual element-based flux recovery on quadtree
%              https://arxiv.org/abs/2006.05585
%   
%
% See also fixorder, fixorder3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

elemVertexNumber = cellfun('length',elem);
maxNv = max(elemVertexNumber);
minNv = min(elemVertexNumber);
for Nv = minNv:maxNv
    isNv = (elemVertexNumber == Nv);
    if any(isNv)
        elemNv = cell2mat(elem(isNv,:));
        % number of elements sharing the same # of vertices
        NelemNv = sum(isNv);
        
        x1 = reshape(node(elemNv,1),NelemNv,Nv);
        y1 = reshape(node(elemNv,2),NelemNv,Nv);
        xc = mean(x1,2)-eps; yc = mean(y1,2); % simple mean with a shift
        nodeNv = node(elemNv',:);
        nodeNv2elem = reshape(nodeNv', [2, Nv, NelemNv]);
        nodeNv2elem = permute(nodeNv2elem, [3 1 2]);
        vert2Centroid = nodeNv2elem - repmat([xc, yc], [1 1 Nv]);
        vert2CentroidAngle = zeros(NelemNv,Nv);
        for j = 1:Nv
            vert2CentroidAngle(:,j) = atan2(vert2Centroid(:,2,j), vert2Centroid(:,1,j));
        end
        [~, ixSorted] = sort(vert2CentroidAngle,2);
        ixSorted = repmat(0:Nv:Nv*(NelemNv-1), [Nv 1]) + ixSorted';
        elemNv = elemNv';
        elemNvSorted = reshape(elemNv(ixSorted(:)), [Nv,NelemNv])';
        elem(isNv) = num2cell(elemNvSorted,2);
    end
end

if nargin > 2
    if strcmp(varargin{1}, 'Clockwise')
        elem = cellfun(@(x) fliplr(x), elem, 'UniformOutput',false);
    end
end
end