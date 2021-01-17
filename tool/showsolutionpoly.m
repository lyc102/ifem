function h = showsolutionpoly(node,elem,u,varargin)
%% SHOWSOLUTIONPOLY plots the solution u on a polygonal mesh in 2-D.
%
%    showmeshpoly(node,elem, u) displays u on a polygonal mesh. elem{j} is a cell 
%    array with j-th cell representing the j-th element.
%
%    showmeshpoly(node,elem,u,viewangle) changes the display angle. The
%    deault view angle on planar meshes is view(2) and view(3) for surface
%    meshes. 
%
%   Example:
%
%     load polygonunstruct
%     u = @(p) sin(2*pi*p(:,1)).*cos(2*pi*p(:,2));
%     uI = u(node);
%     showsolutionpoly(node,elem,uI);
%     
%   See also showmesh, showboundary3, showmeshpoly.
%
% Added by Shuhao Cao.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

V = [node, u];
elemVertexNumber = cellfun('length',elem); % the number of vertices per element
for Nv = min(elemVertexNumber):max(elemVertexNumber)
    idx = (elemVertexNumber == Nv); % index of elements with Nv vertices
    elemNv = cell2mat(elem(idx));
    if size(elemNv,2) == 1
        elemNv = reshape(elemNv,[Nv,nnz(idx)])';
    end
    h = patch('faces',elemNv,'vertices',V);
    hold on;
    set(h,    'FaceVertexCData',u,...
              'Facecolor','interp',...
              'Edgecolor','k',...
              'Edgealpha',1/log(size(node,1)),...
              'FaceAlpha',0.65);
   %%
   if (nargin>3) && ~isempty(varargin) % set display property
       if isnumeric(varargin{1})
           view(varargin{1});
           if nargin>4
               set(h,varargin{2:end});
           end
       else
           set(h,varargin{1:end});
       end
   end
end
view(3); axis tight; grid on;
camproj('perspective');
hold off
%%

