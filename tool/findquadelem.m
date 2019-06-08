function findquadelem(node,elem,range,varargin)
%% FINDQUADELEM highlights some elements
%
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on

if (nargin==2) || isempty(range) || (ischar(range) && strcmp(range,'all'))
    range = (1:size(elem,1))'; 
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
center = ( node(elem(range,1),:) + node(elem(range,2),:) ...
         + node(elem(range,3),:) + node(elem(range,4),:))/4;
if length(range) < size(elem,1)
    x = reshape(node(elem(range,:),1),size(range,1), size(elem,2));
    y = reshape(node(elem(range,:),2),size(range,1), size(elem,2));
    h = patch(x,y,'y');
    if nargin > 3
        if strcmp(varargin{1},'noindex') || strcmp(varargin{1},'index')
           if size(varargin,2)>=2
                set(h,varargin{2:end});
           end
        else
            set(h,varargin{:});
        end
    end
end
if (nargin <=3) || ~(strcmp(varargin{1},'noindex'))
    if size(node,2) == 2
        plot(center(:,1),center(:,2),'o','LineWidth',1,'MarkerEdgeColor','k',...
         'MarkerFaceColor','y','MarkerSize',20);    
        text(center(:,1)-0.015,center(:,2),int2str(range),'FontSize',12,...
        'FontWeight','bold','Color','k');
    elseif size(node,2) == 3 % surface mesh
        plot3(center(:,1),center(:,2),center(:,3),'o','LineWidth',1,'MarkerEdgeColor','k',...
         'MarkerFaceColor','y','MarkerSize',20);    
        text(center(:,1)-0.015,center(:,2),center(:,3),int2str(range),'FontSize',12,...
        'FontWeight','bold','Color','k');        
    end
end
hold off