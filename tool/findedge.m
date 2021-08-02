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
%   See also findelem3, findnode3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on
% set up range
if (nargin==2) || isempty(range) || (ischar(range) && strcmp(range,'all'))
    range = (1:size(edge,1))'; 
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
varidx = 1;
if nargin > 3 && strcmp(varargin{1},'noindex') % not showing index
    varidx = 2;
else
    if ndim == 2
        h = plot(midEdge(:,1),midEdge(:,2),'o','LineWidth',1,'MarkerEdgeColor','k',...
             'MarkerFaceColor',[0.6 0.5 0.8],'MarkerSize',20);
    elseif ndim == 3
        h = plot3(midEdge(:,1),midEdge(:,2),midEdge(:,3),'o','LineWidth',1,'MarkerEdgeColor','k',...
             'MarkerFaceColor',[0.6 0.5 0.8],'MarkerSize',20);
    end    
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
    if (strcmp(varargin{1},'noindex') ...
           && (strcmp(varargin{varidx},'vec') || strcmp(varargin{varidx},'rotvec')...
            || strcmp(varargin{varidx},'draw')) )
        startidx = 3;
    else
        startidx = 2;
    end
    if  length(varargin) >= startidx
        set(h,varargin{startidx:end});
    end
end
if (nargin <= 4) || ~(strcmp(varargin{1},'noindex'))
    if ndim == 2
        text(midEdge(:,1)-0.025,midEdge(:,2),int2str(range), ...
             'FontSize',12,'FontWeight','bold','Color','k');
    elseif ndim == 3
        text(midEdge(:,1)-0.025,midEdge(:,2),midEdge(:,3),int2str(range), ...
         'FontSize',12,'FontWeight','bold','Color','k');
    end
end