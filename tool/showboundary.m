function showboundary(node,elem,expr,varargin)
%% SHOWBOUNDARY boundary mesh plot
%
%    showboundary(node,elem) displays the boundary surface mesh of a
%    3-dimensional tetrahderon mesh given by node and elem matrices. The
%    boundary is found by the function findboundary3.
%
%    showboundary(node,elem,expr) displays the boundary of parts of the
%    mesh specificed by the expression. For example,
%    showboundary3(node,elem,'~(x>=0 & y>=0 & z>=0)') only shows the boundary
%    mesh for tetrahedrons not in the first quadrant. 
%
%    showboundary(node,elem,expr,viewangle) changes the display angle. The
%    deault view angle is view(3). 
%
%    showboundary(node,elem,expr,'param','value','param','value'...)
%    allows additional patch param/value pairs to be used when displaying
%    the mesh. For example, the default transparency parameter is set to
%    0.5. You can overwrite this value by using the param pair
%    ('FaceAlpha', value). The value has to be a number between 0 and 1.
%    Other parameters include: 'Facecolor', 'Edgecolor' etc.
%   
%    To display the tetrahedron mesh in 3-d, use showmesh3. Notice that the
%    3-d graphical visulization is slow for large mesh data. 
%
%   Example:
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     [node,elem] = uniformbisect3(node,elem);
%     subplot(1,2,1);
%     showboundary3(node,elem);
%     subplot(1,2,2);
%     showboundary3(node,elem,'~(x>=0 & y>=0 & z>=0)',[59,20]); 
%
%   See also showboundary, showsolution3, showmesh.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (nargin >= 3) && (any(expr))
    x = node(:,1);  y = node(:,2);  z = node(:,3); %#ok<*NASGU>
    incl = find(eval(expr));
    elem = elem(any(ismember(elem,incl),2),:);
end
[bdNode, bdFace] = findboundary3(elem); %#ok<*ASGLU>
if isempty(varargin)
    showmesh(node,bdFace);
else
    showmesh(node,bdFace,varargin{1:end});
end

