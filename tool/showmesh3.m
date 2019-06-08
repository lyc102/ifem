function showmesh3(node,elem,expr,varargin)
%% SHOWMESH3 displays a tetrahedron mesh in 3-D.
%
%    showmesh3(node,elem) displays a 3-dimensional tetrahderon mesh given
%    by node and elem matrices; see <a href="matlab:ifem('meshdoc')">meshdoc</a> for the data structure:
%    node and elem.
%
%    showmesh3(node,elem,expr) displays parts of the mesh specificed by the
%    expression. For example, showmesh3(node,elem,'~(x>=0 & y>=0)') only
%    shows the tetrahedron not in the first quadrant. 
%
%    showmesh3(node,elem,expr,'param','value','param','value'...) allows
%    additional patch param/value pairs to be used when displaying the
%    mesh. For example, the default transparency parameter is set to 0.5.
%    You can overwrite this value by using the param pair ('FaceAlpha',
%    value). The value has to be a number between 0 and 1. Other parameters
%    include: 'Facecolor', 'Edgecolor' etc.
%   
%    For meshes with large data, the 3-D graphics is very slow. You may use
%    <a href="matlab:help showboundary3">showboundary3</a> to display the boundary surface mesh only.
%
%   Example:
%     % A mesh for a cube
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     [node,elem] = uniformbisect3(node,elem);
%     subplot(1,2,1);
%     showmesh3(node,elem); pause(1)
%     subplot(1,2,2);
%     showmesh3(node,elem,'~(x>=0 & y>=0 & z>=0)','FaceAlpha',0.25); 
%     axis on; view([59,20])
%
%   See also showboundary3, showsolution3, showmesh.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (nargin >= 3) && (any(expr))
    x = node(:,1);  y = node(:,2);  z = node(:,3); %#ok<NASGU>
    incl = find(eval(expr));
    elem = elem(any(ismember(elem,incl),2),:);
end
h = tetramesh(elem(:,1:4),node,ones(size(elem,1),1));
set(h,'facecolor',[0.35 0.75 0.35],'edgecolor','k');
if nargin > 3 
    set(h,varargin{1:end})
else % default display properties
    set(h,'FaceAlpha',0.4);
end
view(3);
axis off; axis equal; axis tight