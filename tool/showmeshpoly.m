function h = showmeshpoly(node,elem,varargin)
%% showmeshpoly displays a polygonal/polyhedron mesh in 2-D or 3-D.
%
%    showmeshpoly(node,elem) displays a polygonal mesh. elem{j} is a cell 
%    array with j-th cell representing the j-th element.
%
%    showmeshpoly(node,elem,viewangle) changes the display angle. The
%    deault view angle on planar meshes is view(2) and view(3) for surface
%    meshes. 
%        
%    showmeshpoly(node,elem,'param','value','param','value'...) allows
%    additional patch param/value pairs to be used when displaying the
%    mesh.  For example, the default transparency parameter for a surface
%    mesh is set to 0.75. You can overwrite this value by using the param
%    pair ('FaceAlpha', value). The value has to be a number between 0 and
%    1. Other parameters include: 'Facecolor', 'Edgecolor' etc. These
%    parameters are mostly used when displaying a surface mesh.
%
%   Example:
%
%     load polygonunstruct
%     showmeshpoly(node,elem);
%     
%   See also showmesh, showboundary3, showsolutionvem.
%
% S. Cao modified from Long Chen's iFEM vem codes. 2016

dim = size(node,2);

if (dim==2) % 2D mesh
    elemVertexNumber = cellfun('length',elem);% the number of vertices per element
    for Nv = min(elemVertexNumber):max(elemVertexNumber)
        idx = (elemVertexNumber == Nv); % index of elements with Nv vertices
        elemNv = cell2mat(elem(idx));
        if size(elemNv,2)==1
            elemNv = reshape(elemNv,Nv,nnz(idx))';
        end
        h = patch('faces',elemNv,'vertices',node);
        hold on;
        set(h,'FaceColor',[0.5, 0.8, 0.8],...
              'Edgecolor','k',...
              'Edgealpha',1/nthroot(size(node,1),3),...
              'FaceAlpha',0.75);
    end
    view(2); axis equal off tight;
end

%% needs modfication
if (dim==3) % 3D mesh
    elemVertexNumber = cellfun('length',elem);% the number of vertices per element
    %to-do
    if min(elemVertexNumber)>=4 % vertex # in any elem >=4
        idx4 = (elemVertexNumber == 4); % index of elements with 4 vertices
        elemNv4 = cell2mat(elem(idx4));
        if size(elemNv,2)==1
            elemNv4 = reshape(elemNv4,4,nnz(idx4))';
        end
        normalVec1 = mycross(node(elemNv4(:,1),:) - node(elemNv4(:,3),:),...
            node(elemNv4(:,1),:) - node(elemNv4(:,2),:));
        normalVec2 = mycross(node(elemNv4(:,1),:) - node(elemNv4(:,3),:),...
            node(elemNv4(:,1),:) - node(elemNv4(:,4),:));
        isCoplanar = mycross(normalVec1,normalVec2);
        if norm(isCoplanar)<1e-10 % surface mesh
            showmeshpoly(node,elem,varargin{:});
        end
    elseif max(elemVertexNumber)<=3  % surface triangles
        showmesh(node,elem,varargin{:})
    else 
        showmesh3(node,elem,varargin{:});
    end
end

%%
if (nargin>2) && ~isempty(varargin) % set display property
    if isnumeric(varargin{1})
        view(varargin{1});
        if nargin>3
            set(h,varargin{2:end});
        end
    else
        set(h,varargin{1:end});        
    end
end
