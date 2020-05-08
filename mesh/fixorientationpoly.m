function elem = fixorientationpoly(node,elem,varargin)
%% FIXORIENTATIONPOLY set all polygon clockwise or counter-clockwise
% 
%   elem = fixorientationpoly(node,elem) reorder the elem to orient
%   clockwise similar to the quad element set up. 
%   elem{:} has to be a cell array.
%
%
% See also fixorder, fixorder3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

elem = cellfun(@(x) reorder(x), elem,'UniformOutput',false);
if nargin > 2
    if strcmp(varargin{1}, 'CounterClockwise')
        elem = cellfun(@(x) fliplr(x), elem, 'UniformOutput',false);
    end
end

%%
    function idx = reorder(idx)
        [xc, yc] = centroid(node(idx,:));
        vecx = node(idx,:) - [xc, yc];
        [~, idxOrdered] = sort(atan2(vecx(:,1), vecx(:,2)));
        idx = idx(idxOrdered);
    end

    function [xc, yc] = centroid(p)
        x1 = p(:,1);
        y1 = p(:,2);
        x2 = circshift(x1,-1);
        y2 = circshift(y1,-1);
        bdIntegral = x1.*y2 - y1.*x2;
        area = sum(bdIntegral)/2;
        xc = sum((x1+x2).*bdIntegral) /6/area;
        yc = sum((y1+y2).*bdIntegral) /6/area;
    end
end