function showresult3(node,elem,u,expr,varargin)
%% SHOWRESULT3 display the mesh and the solution in 3-D
%
%  showresult3(node,elem,u,viewangle) displays the mesh and the solution in
%  one figure. The left one is the mesh, the middle
%  one is the contour of the solution, and the right one is the graph of
%  the function. The last viewangle is used to adjust the view angle of the
%  graph of the function.
%
%  Example:
%
% See also showrate, showmesh, showsolution
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('expr','var'), expr = []; end
set(gcf,'Units','normal'); 
set(gcf,'Position',[0.25,0.25,0.6,0.4]);
% show mesh
if size(elem,1) < 1e6
    subplot(1,2,1); 
    showboundary3(node,elem,expr,varargin{:});
else
    subplot(1,2,1);
    title('The mesh is too dense to display')
end
% show solution
subplot(1,2,2);
if size(u,1) == size(node,1)
    showsolution3(node,elem,u,expr,varargin{:});
    colorbar;
end
pause(0.05)