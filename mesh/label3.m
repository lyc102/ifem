function [elem,bdFlag] = label3(node,elem,t,bdFlag)
%% LABEL3 label the refinement edge of a 3-D triangulation
%
% elem = label3(node,elem) return a triangulation such that elem(t,[1 2]) is
% the longest edge of t.
% 
% [elem,bdFlag] = label3(node,elem,bdFlag) also switch the boundary edges
% accordingly.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (nargin==2)
    t=1:size(elem,1); 
    bdFlag = [];
end

%% Compute edge lengths
allEdge = [elem(t,[1 2]); elem(t,[1 3]); elem(t,[1 4]); ...
           elem(t,[2 3]); elem(t,[2 4]); elem(t,[3 4])];
allEdgeLength = sum((node(allEdge(:,1),:)-node(allEdge(:,2),:)).^2,2);
allEdgeLength = (1 + 0.1*rand(size(allEdge,1),1)).*allEdgeLength;
elemEdgeLength = reshape(allEdgeLength,length(t),6);
[tempvar,idx] = max(elemEdgeLength,[],2);

%% Reorder the vertices of elem
elem(t((idx==2)),1:4) = elem(t((idx==2)),[3 1 2 4]);
elem(t((idx==3)),1:4) = elem(t((idx==3)),[1 4 2 3]);
elem(t((idx==4)),1:4) = elem(t((idx==4)),[2 3 1 4]);
elem(t((idx==5)),1:4) = elem(t((idx==5)),[2 4 3 1]);
elem(t((idx==6)),1:4) = elem(t((idx==6)),[4 3 2 1]);

%% Reorder the boundary flags
if (nargin==4) && (~isempty(bdFlag))
    bdFlag(t((idx==2)),1:4) = bdFlag(t((idx==2)),[3 1 2 4]);
    bdFlag(t((idx==3)),1:4) = bdFlag(t((idx==3)),[1 4 2 3]);
    bdFlag(t((idx==4)),1:4) = bdFlag(t((idx==4)),[2 3 1 4]);
    bdFlag(t((idx==5)),1:4) = bdFlag(t((idx==5)),[2 4 3 1]);
    bdFlag(t((idx==6)),1:4) = bdFlag(t((idx==6)),[4 3 2 1]);
else
    bdFlag = [];
end