function h = showsolutionpoly(node,elem,u,varargin)
%% SHOWSOLUTIONPOLY plots the solution u on a polygonal mesh in 2-D.
%
%    showmeshpoly(node,elem, u) displays u on a polygonal mesh. elem{j} is a cell 
%    array with j-th cell representing the j-th element.
%    
%    u can be either:
%    - Nodal values: size(u,1) = size(node,1) for continuous interpolation
%    - Element values: size(u,1) = length(elem) for piecewise constant
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

if (size(u,1) == length(elem))
    % For piecewise constant: use flat shading with element-based colors
    V = node;
    elemVertexNumber = cellfun('length',elem);
    elementIdx = 1;
    
    for Nv = min(elemVertexNumber):max(elemVertexNumber)
        idx = (elemVertexNumber == Nv);
        elemNv = cell2mat(elem(idx));
        if size(elemNv,2) == 1
            elemNv = reshape(elemNv,[Nv,nnz(idx)])';
        end
        
        % Get element values for current group
        uElem = u(idx);
        
        h = patch('faces',elemNv,'vertices',V);
        hold on;
        set(h,    'FaceVertexCData',uElem,...
                  'Facecolor','flat',...
                  'Edgecolor','k',...
                  'Edgealpha',1/log(size(node,1)),...
                  'FaceAlpha',0.65);
        elementIdx = elementIdx + nnz(idx);
    end
elseif (size(u,1) == size(node,1))
    % Original code for nodal values
    V = [node, u];
    elemVertexNumber = cellfun('length',elem);
    
    for Nv = min(elemVertexNumber):max(elemVertexNumber)
        idx = (elemVertexNumber == Nv);
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
    end
else
    error('Invalid size of u. It must match either number of nodes (%d) or number of elements (%d)', ...
          size(node,1), length(elem));
end

%% Handle optional arguments
if (nargin>3) && ~isempty(varargin)
    if isnumeric(varargin{1})
        view(varargin{1});
        if nargin>4
            set(h,varargin{2:end});
        end
    else
        set(h,varargin{1:end});
    end
end

view(3); axis tight; grid on;
camproj('perspective');
hold off

