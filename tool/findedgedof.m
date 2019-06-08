function findedgedof(node,edge,range,varargin)
%% FINDEDGEDOF highlights some dofs associated to edges
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on
N = size(node,1);
% set up range
if (nargin<=2) || (isempty(range)) 
    range = (1:size(edge,1))'+ N; 
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
dofrange = range;
range = range - N;
midEdge = (node(edge(range,1),:)+node(edge(range,2),:))/2;
if length(range) < size(edge,1)
    h = line([node(edge(range,1),1); node(edge(range,2),1)],...
         [node(edge(range,1),2); node(edge(range,2),2)],...
         'LineWidth',1.5,'Color','r');
else
    h = plot(midEdge(:,1),midEdge(:,2),'r.','MarkerSize', 22);
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
    text(midEdge(:,1)+0.025,midEdge(:,2)-0.025,int2str(dofrange), ...
         'FontSize',12,'FontWeight','bold','Color','k');
end
% TODO: help and example