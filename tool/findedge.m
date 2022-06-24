function findedge(node,edge,range,varargin)
%% FINDEDGE highlights edges
%
%    FINDEDGE(node,edge,range) finds all elements whose indices are in the
%    range array by displaying these elements in yellow.
%
%    FINDEDGE(node,edge) finds all elements.
%   
%    FINDEDGE(node,edge,1,'draw') skip the display of indices.
%
%    FINDEDGE(node,edge,range,'param','value','param','value'...) allows
%    additional patch param/value pairs to highlight the elements.
%    
% Example:
%     [node,elem] = squaremesh([0,1,0,1],1/2);
%     subplot(1,2,1);
%     showmesh(node,elem);
%     T = auxstructure(elem);
%     findedge(node,T.edge,5,'noindex','MarkerFaceColor','r','MarkerSize',24);
%     subplot(1,2,2);
%     showmesh(node,elem);
%     findedge(node,T.edge);
%
% Example: local labeling of edges
%
%     node = [1,0; 1,1; 0,0];
%     elem = [1 2 3];
%     edge = [2 3; 1 3; 1 2];
%     showmesh(node,elem);
%     findnode(node);
%     findedge(node,edge,'all','MarkerSize',28);
%
%   See also findelem3, findnode3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on
dotColor = 'k.';
% set up range
if (nargin==2) || isempty(range) || (ischar(range) && strcmp(range,'all'))
    range = (1:size(edge,1))'; 
    dotColor = 'r.';
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
% draw edges in red
ndim = size(node,2);
midEdge = (node(edge(range,1),:)+node(edge(range,2),:))/2;
if ndim == 2
    h = plot(midEdge(:,1),midEdge(:,2),dotColor,'MarkerSize',20);
elseif ndim == 3
    h = plot3(midEdge(:,1),midEdge(:,2),midEdge(:,3),dotColor,'MarkerSize',20);
end    
varidx = 1;
if nargin > 3 && strcmp(varargin{1},'noindex') % not showing index
    varidx = 2;
end
if nargin > 3
    if strcmp(varargin{varidx},'draw') 
        switch ndim
            case 2
            h = line([node(edge(range,1),1)'; node(edge(range,2),1)'],...
                     [node(edge(range,1),2)'; node(edge(range,2),2)'],...
                    'LineWidth',2,'Color','r');
            case 3
            h = line([node(edge(range,1),1)'; node(edge(range,2),1)'],...
                     [node(edge(range,1),2)'; node(edge(range,2),2)'],...
                     [node(edge(range,1),3)'; node(edge(range,2),3)'],...
                    'LineWidth',2,'Color','r');
        end
    end
    if strcmp(varargin{varidx},'vec') % plot edge vector
        edgeVec = node(edge(range,2),:) - node(edge(range,1),:);    
        if ndim == 2
            h = quiver(midEdge(:,1),midEdge(:,2),edgeVec(:,1),edgeVec(:,2));
        elseif ndim ==3
            h = quiver3(midEdge(:,1),midEdge(:,2),midEdge(:,3), ...
                        edgeVec(:,1),edgeVec(:,2),edgeVec(:,3));
        end
        set(h,'Linewidth',3)
    end
    if strcmp(varargin{varidx},'rotvec') && ndim == 2 % plot edge normal vector
        edgeVec = node(edge(range,2),:) - node(edge(range,1),:);    
        h = quiver(midEdge(:,1),midEdge(:,2),-edgeVec(:,2),edgeVec(:,1));
        set(h,'Linewidth',3)
    end
    startidx = 1;
    if strcmp(varargin{1},'noindex') 
        startidx = startidx + 1;
    end
    if (strcmp(varargin{varidx},'vec') || strcmp(varargin{varidx},'rotvec')...
            || strcmp(varargin{varidx},'draw'))
        startidx = startidx + 1;
    end
    if  length(varargin) >= startidx
        set(h,varargin{startidx:end});
    end
end
if (nargin <= 4) || ~(strcmp(varargin{1},'noindex'))
    if ndim == 2
        text(midEdge(:,1)+0.025,midEdge(:,2)+0.015,int2str(range), ...
             'FontSize',12,'FontWeight','bold','Color','k');
    elseif ndim == 3
        text(midEdge(:,1)+0.025,midEdge(:,2)+0.015,midEdge(:,3),int2str(range), ...
         'FontSize',12,'FontWeight','bold','Color','k');
    end
end