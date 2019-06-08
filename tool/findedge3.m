function findedge3(node,edge,range,varargin)
%% FINDEDGE3 highlights edges in certain range in 3-D
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on
% set up range
if (nargin<=2) || (isempty(range)) 
    range = (1:size(edge,1))'; 
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
midEdge = (node(edge(range,1),:)+node(edge(range,2),:))/2;
if length(range)<size(edge,1)
    h = plot3([node(edge(range,1),1)'; node(edge(range,2),1)'],...
              [node(edge(range,1),2)'; node(edge(range,2),2)'],...
              [node(edge(range,1),3)'; node(edge(range,2),3)'],...
              'LineWidth',1.5,'Color','r');
else
    h = plot3(midEdge(:,1),midEdge(:,2),midEdge(:,3),'r.','MarkerSize', 18);
end
if nargin > 3
    if strcmp(varargin{1},'noindex') || strcmp(varargin{1},'index')
        startidx = 2;
    else
        startidx = 1;
    end
    set(h,varargin{startidx:end});
end
if (nargin <=3) || ~(strcmp(varargin{1},'noindex'))
    text(midEdge(:,1)+0.025,midEdge(:,2)-0.025,midEdge(:,3),int2str(range), ...
         'FontSize',12,'FontWeight','bold','Color','k');
end
%% TODO: help and example