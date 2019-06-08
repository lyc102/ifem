function showgraph(node,edge)
%% SHOWGRAPH displays a planar graph
%
%    showgraph(node,edge) displays a planar undirected graph.
%
%   See also showmesh, findedge
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

dim = size(node,2);
if dim == 2 % 2-D
    line([node(edge(:,1),1)'; node(edge(:,2),1)'],...
         [node(edge(:,1),2)'; node(edge(:,2),2)'],...
         'LineWidth',1,'Color',[0.125 0.5 0.125]);
    hold on
    plot(node(:,1),node(:,2),'k.', 'MarkerSize', 12);
elseif dim == 3 % 3-D
    line([node(edge(:,1),1)'; node(edge(:,2),1)'],...
         [node(edge(:,1),2)'; node(edge(:,2),2)'],...
         [node(edge(:,1),3)'; node(edge(:,2),3)'],...
         'LineWidth',1,'Color',[0.125 0.5 0.125]);
    hold on
    plot3(node(:,1),node(:,2),node(:,3),'k.', 'MarkerSize', 22);    
end
axis equal; axis off