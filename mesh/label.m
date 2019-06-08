function [elem,bdFlag] = label(node,elem,bdFlag)
%% LABEL label the refinement edge of a 2-D triangulation
%
% elem = label(node,elem) return a tirangulation such that elem(t,1) is
% opposite to the longest edge of t.
%
% [elem,bdFlag] = label(node,elem,bdFlag) also switch the boundary edges
% accordingly.
%  
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Compute length of each edge
totalEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
ve = node(totalEdge(:,1),:) - node(totalEdge(:,2),:);
edgeLength = reshape(sum(ve.^2,2),size(elem,1),3);

%% Switch indices according the edge length
[tempvar,I] = max(edgeLength,[],2); %#ok<*ASGLU>
elem((I==2),[1 2 3]) = elem((I==2), [2 3 1]);
elem((I==3),[1 2 3]) = elem((I==3), [3 1 2]);

%% Switch the boundary edges
if nargin==3
	bdFlag((I==2),[1 2 3]) = bdFlag((I==2), [2 3 1]);
	bdFlag((I==3),[1 2 3]) = bdFlag((I==3), [3 1 2]);
else
	bdFlag = [];
end