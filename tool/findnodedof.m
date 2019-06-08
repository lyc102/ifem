function findnodedof(node,elem,elem2dof)
%% FINDNODEDOF highlights dof associated to nodes.
%
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

hold on
range = 1:size(node,1);
plot(node(range,1),node(range,2),'k.', 'MarkerSize', 18);
text(node(elem(:),1)+0.025,node(elem(:),2)+0.025,int2str(elem2dof(:)), ...
     'FontSize',16,'FontWeight','bold');
hold off