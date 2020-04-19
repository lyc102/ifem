function [elem,bdFlag,HB] = label3(node,elem,markedElem,bdFlag)
%% LABEL3 label the refinement edge of a 3-D triangulation
%
% elem = label3(node,elem) return a triangulation such that elem(t,[1 2]) is
% the longest edge of t.
% 
% [elem,bdFlag,HB] = label3(node,elem,markedElem,bdFlag) also switch the
% boundary flags but for marked elements only. HB is needed when use
% coarsen3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (nargin==2)
    markedElem = 1:size(elem,1); 
    bdFlag = [];
end
if strcmp(markedElem,'all'), markedElem = (1:size(elem,1))'; end

%% Compute edge lengths
allEdge = [elem(markedElem,[1 2]); elem(markedElem,[1 3]); elem(markedElem,[1 4]); ...
           elem(markedElem,[2 3]); elem(markedElem,[2 4]); elem(markedElem,[3 4])];
allEdgeLength = sum((node(allEdge(:,1),:)-node(allEdge(:,2),:)).^2,2);
allEdgeLength = (1 + 0.1*rand(size(allEdge,1),1)).*allEdgeLength;
elemEdgeLength = reshape(allEdgeLength,length(markedElem),6);
[tempvar,idx] = max(elemEdgeLength,[],2);

%% Reorder the vertices of elem
elem(markedElem((idx==2)),1:4) = elem(markedElem((idx==2)),[3 1 2 4]);
elem(markedElem((idx==3)),1:4) = elem(markedElem((idx==3)),[1 4 2 3]);
elem(markedElem((idx==4)),1:4) = elem(markedElem((idx==4)),[2 3 1 4]);
elem(markedElem((idx==5)),1:4) = elem(markedElem((idx==5)),[2 4 3 1]);
elem(markedElem((idx==6)),1:4) = elem(markedElem((idx==6)),[4 3 2 1]);

%% Reorder the boundary flags
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlag(markedElem((idx==2)),1:4) = bdFlag(markedElem((idx==2)),[3 1 2 4]);
    bdFlag(markedElem((idx==3)),1:4) = bdFlag(markedElem((idx==3)),[1 4 2 3]);
    bdFlag(markedElem((idx==4)),1:4) = bdFlag(markedElem((idx==4)),[2 3 1 4]);
    bdFlag(markedElem((idx==5)),1:4) = bdFlag(markedElem((idx==5)),[2 4 3 1]);
    bdFlag(markedElem((idx==6)),1:4) = bdFlag(markedElem((idx==6)),[4 3 2 1]);
else
    bdFlag = [];
end

%% Set up HB for bisection and coarsening
N = size(node,1);
HB = zeros(N,4);
HB(1:N,1:3) = repmat((1:N)',1,3);