function findelem3(node,elem,range,varargin)
%% FINDELEM3 highlights the elements in the range
%
%    findelem3(node,elem,2)
%    findelem3(node,elem,[0 1 0])
%    findelem3(node,elem)
%
% Input: 
%    range is a vector containing index of element. It could be also logic.
%    If range is skipped, findelem will find all elements.
%
%   See also findelem, findnode3, findedge.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on
if (nargin==2)
    range = (1:size(elem,1))'; 
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
dimension = size(elem,2)-1;
if dimension == 3
    center = (node(elem(range,1),:) + node(elem(range,2),:) + ...
              node(elem(range,3),:) + node(elem(range,4),:))/4;
else
    center = (node(elem(range,1),:) + node(elem(range,2),:) + ...
              node(elem(range,3),:))/3;
end    
plot3(center(:,1),center(:,2),center(:,3),'o','LineWidth',1,...
                          'MarkerEdgeColor','k','MarkerFaceColor','y',...
                          'MarkerSize',20)
%                            'MarkerFaceColor','w',...
text(center(:,1)-0.02,center(:,2)-0.02,center(:,3),...
         int2str(range),'FontSize',14,'FontWeight','bold');
hold off;