function showresult(node,elem,u,viewangle)
%% SHOWRESULT display the mesh and the solution 
%
%  showresult(node,elem,u,viewangle) displays the mesh and the solution in
%  one figure. The left one is the mesh, the middle
%  one is the contour of the solution, and the right one is the graph of
%  the function. The last viewangle is used to adjust the view angle of the
%  graph of the function.
%
%  Example:
%     f = inline('sin(2*pi*x).*cos(2*pi*y)');
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     for k = 1:4
%         [node,elem] = uniformrefine(node,elem);
%     end
%     u = f(node(:,1),node(:,2));
%     showresult(node,elem,u,[-62,58]);
%
% See also showrate, showmesh, showsolution
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Plot mesh
set(gcf,'Units','normal'); 
set(gcf,'Position',[0.25,0.25,0.7,0.25]);
if size(elem,1) < 6e3
    subplot(1,3,1); 
    showmesh(node,elem); 
%     pause(0.05)
else
    subplot(1,3,1);
    title('The mesh is too dense to display')
end
%% Plot solution on plane
subplot(1,3,2); 
showsolution(node,elem,u,2);
colorbar;

%% Plot graph of the solution
subplot(1,3,3); 
showsolution(node,elem,u);
if nargin>3
    view(viewangle);
end
pause(0.05)