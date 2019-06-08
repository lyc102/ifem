function uI = Lagrangeinterpolate(u,node,elem,elemType,varargin)
%% LAGRANGEINTERPOLATE interpolate to Lagrange elements

if ~exist('elemType','var')
    elemType = 'P1';
end
dim = size(elem,2) - 1;
switch elemType
    case 'P0'     % piecewise constant function P1 element
        if dim == 2
            barycenter = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
        elseif dim == 3
            barycenter = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:) +node(elem(:,4),:))/4;
        end
        uI = u(barycenter);
    case 'P1'     % piecewise linear function P1 element
        uI = u(node); % nodal interpolation
    case 'CR'     % piecewise linear function CR element
        if dim == 2
            edge = varargin{1};
            uI = u((node(edge(:,1),:)+node(edge(:,2),:))/2);                
        elseif dim == 3
            face = varargin{1};
            uI = u((node(face(:,1),:)+node(face(:,2),:)+node(face(:,3),:))/3);                
        end
    case 'P2'     % piecewise quadratic function
        edge = varargin{1};
        uI = u([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
    case 'P3' 
        edge = varargin{1};
        Ne = size(edge,1);
        N  = size(node,1);
        NT = size(elem,1);
        uI = u(node);
        uI(N+1:2:N+2*Ne-1) = u(node(edge(:,1),:)+1/3*(node(edge(:,2),:)-node(edge(:,1),:)));          
        uI(N+2:2:N+2*Ne) = u(node(edge(:,1),:)+2/3*(node(edge(:,2),:)-node(edge(:,1),:)));
        uI(N+2*Ne+1:N+2*Ne+NT) = u((node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3);
    case 'WG'     % weak Galerkin element
        if dim == 2
            edge = varargin{1};
            center = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
            midedge = (node(edge(:,1),:)+node(edge(:,2),:))/2;
            uI = u([center; midedge]);
        elseif dim == 3
            face = varargin{1};
            center = (node(elem(:,1),:)+node(elem(:,2),:)...
                     +node(elem(:,3),:)+node(elem(:,4),:))/4;
            midface = (node(face(:,1),:)+node(face(:,2),:)+node(face(:,3),:))/3;
            uI = u([center; midface]);
        end            
end
    
