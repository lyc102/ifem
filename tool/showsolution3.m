function showsolution3(node,elem,u,expr,varargin)
%% SHOWSOLUTION3 plots the solution u on a tetrahedron mesh in 3-D.
%
%    showsolution3(node,elem,u) displays the functoin u on a topological
%    2-dimensional mesh given by node and elem matrices. The function u
%    could be piecewise constant or piecewise linear. 
%
%    showsolution3(node,elem,u,expr) displays the function u on the
%    boundary of parts of the mesh specificed by the expression. For
%    example, showsoluiton3(node,elem,'z==0') will show the function on the
%    cross section of z=0. 
%
%    showsolution3(node,elem,u,viewangle) changes the display angle. The
%    deault view angle on planar meshes is view(2) and view(3) for surface
%    meshes. 
%
%    showsolution3(node,elem,u,expr,'param','value','param','value'...) allows
%    additional patch param/value pairs to be used when displaying the
%    mesh. 
%
%   Example:
%     f = inline('x.^2 + y.^2 + z.^2');
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     for k=1:4
%       [node,elem] = uniformbisect3(node,elem);
%     end
%     u = f(node(:,1),node(:,2),node(:,3));
%     subplot(1,2,1);
%     showsolution3(node,elem,u);
%     subplot(1,2,2);
%     showsolution3(node,elem,u,'~(x>0 & y>0)',[139,16],'EdgeColor','k');
%
%   See also showmesh, showsolution3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (nargin >= 4) && (any(expr))
    x = node(:,1);  y = node(:,2);  z = node(:,3); %#ok<*NASGU>
    incl = find(eval(expr));
    elem = elem(any(ismember(elem,incl),2),:);
end
[bdNode, bdFace] = findboundary3(elem); %#ok<ASGLU>
if isempty(varargin)
    showsolution(node,bdFace,u,'EdgeColor','k');
elseif (nargin == 5) && isnumeric(varargin{1})
    showsolution(node,bdFace,u,'EdgeColor','k');
    view(varargin{1});
else
    showsolution(node,bdFace,u,varargin{1:end});
end