function h = showvector3(p,u,option,varargin)
%% showvector3 plot the 3D vector fields using coneplot or quiver3
%
%  showvector3(P,U) displays vector field U associated with points P
%
%  showvector3(P,U,option) displays vector field U associated with points P
%  with various option, see below
%
%  showvector3(P,U,option,'PropertyName',PropertyValue,...) specifies 
%  property name and property value pairs for the quivergroup/patch
%  objects the function creates
% 
% p: Points information, (Number of Points)*3 array
% u: Vector at each point, (Number of Points)*3 array
%
% option.stepSize: the space between each arrows
% option.layer: the number of layers in the Z-direction
% option.scale: scale of the cones or arrows
% option.interp: 'linear' or 'nearest'
% option.color: color can be the magnitude of u
% option.plot: 'coneplot' or 'quiver3'
%
% varargin: see MATLAB help for LineSpec/Figure/quivergroup Properties
% 
% See also quiver3, coneplot, trisurf, showvector3Demo.
% S. Cao 2012

%% Options
if ~exist('option','var')
    option = struct;
end

%% Bounding box for vectors
xmin = min(p(:,1));
xmax = max(p(:,1));
ymin = min(p(:,2));
ymax = max(p(:,2));
zmax = max(p(:,3));
zmin = min(p(:,3));
allmin = min([xmin;ymin;zmin]);
allmax = max([xmax;ymax;zmax]);

if isfield(option,'stepSize')  && isnumeric(option.stepSize)
    step = option.stepSize;
else
    step = (allmax - allmin)/15;
end
[x,y,z] = meshgrid(allmin:step:allmax);

%% Number of layer representation
if isfield(option,'layer') && isnumeric(option.layer)
    if option.layer == 1
        [sx,sy,sz] = meshgrid(xmin:2*step:xmax,...
            ymin:2*step:ymax,(zmin+zmax)/2);
    else
        layerStep = (zmax - zmin)/(option.layer-1);
        [sx,sy,sz] = meshgrid(xmin:2*step:xmax,...
            ymin:2*step:ymax,zmin:layerStep:zmax);
    end
else
    option.layer = 4;
    layerStep = (zmax - zmin)/(option.layer-1);
    [sx,sy,sz] = meshgrid(xmin:2*step:xmax,...
        ymin:2*step:ymax,zmin:layerStep:zmax);
end

%% Interpolation of the data onto regular grid
ux = u(:,1);
uy = u(:,2);
uz = u(:,3);

xc = p(:,1);
yc = p(:,2);
zc = p(:,3);

if ~isfield(option,'interp') || ~ischar(option.interp)
    option.interp = 'linear';
end

if ~isfield(option,'plot') || ~ischar(option.plot)
    option.plot = 'coneplot';
end

switch lower(option.plot(1))
    case 'c'
        switch lower(option.interp(1))
            case 'l'
                uxinterp = gridinterp3(xc,yc,zc,ux,x,y,z,'linear');
                uyinterp = gridinterp3(xc,yc,zc,uy,x,y,z,'linear');
                uzinterp = gridinterp3(xc,yc,zc,uz,x,y,z,'linear');
            case 'n'
                uxinterp = gridinterp3(xc,yc,zc,ux,x,y,z,'nearest');
                uyinterp = gridinterp3(xc,yc,zc,uy,x,y,z,'nearest');
                uzinterp = gridinterp3(xc,yc,zc,uz,x,y,z,'nearest');
        end
        
    case 'q'
        switch lower(option.interp(1))
            case 'l'
                uxinterp = gridinterp3(xc,yc,zc,ux,sx,sy,sz,'linear');
                uyinterp = gridinterp3(xc,yc,zc,uy,sx,sy,sz,'linear');
                uzinterp = gridinterp3(xc,yc,zc,uz,sx,sy,sz,'linear');
            case 'n'
                uxinterp = gridinterp3(xc,yc,zc,ux,sx,sy,sz,'nearest');
                uyinterp = gridinterp3(xc,yc,zc,uy,sx,sy,sz,'nearest');
                uzinterp = gridinterp3(xc,yc,zc,uz,sx,sy,sz,'nearest');
        end
        
    otherwise
end

%% Show the vector field

if isfield(option,'color')
    if size(option.color) == size(uxinterp)
        color = option.color;
    end
else
    % default color uses magnitude of u
    magnitude = sqrt(uxinterp.^2 + uyinterp.^2 + uzinterp.^2);
    color = magnitude;
end


switch lower(option.plot(1))
    case 'c'
        if isfield(option,'scale')
            scale = option.scale;
            if isnumeric(scale) && (scale >= 5)
                scale = 5;
            end
            h = coneplot(x,y,z,uxinterp,uyinterp,uzinterp,...
                sx(:),sy(:),sz(:),scale,color);
        else
            h = coneplot(x,y,z,uxinterp,uyinterp,uzinterp,...
                sx(:),sy(:),sz(:),color);
        end
        set(h,'EdgeColor','none','faceAlpha',0.6);
    case 'q'
        if isfield(option,'scale')
            scale = option.scale;
            if isnumeric(scale) && (scale >= 5)
                scale = 5;
            end
            h = quiver3(sx,sy,sz,uxinterp,uyinterp,uzinterp,scale);
        else
            h = quiver3(sx,sy,sz,uxinterp,uyinterp,uzinterp);
        end
        set(h,'LineWidth',2,'MaxHeadSize',0.5);        
end

%% Setting the view
axis([xmin xmax ymin ymax zmin zmax]); 
if (option.layer == 1)
    view(2);
else
    view(3); grid on; camproj('perspective');
end

if nargin>3
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


function w = gridinterp3(x,y,z,v,xi,yi,zi,method)
%% gridinterp3 Data gridding and hyper-surface fitting for 3-dimensional data.
% Simplified version of TriScatteredInterp

x = x(:); y=y(:); z=z(:); v = v(:);
X = [x y z];

% Sort (x,y,z) so duplicate points can be averaged before passing to delaunay

[X, ind] = sortrows(X);
v = v(ind);
ind = all(diff(X)'==0);
if any(ind)
  ind = [0 ind];
  ind1 = diff(ind);
  fs = find(ind1==1);
  fe = find(ind1==-1);
  if fs(end) == length(ind1) % add an extra term if the last one start at end
     fe = [fe fs(end)+1];
  end
  
  for i = 1:length(fs)
    % averaging v values
    v(fe(i)) = mean(v(fs(i):fe(i)));
  end
  X = X(~ind(2:end),:);
  v = v(~ind(2:end));
end

switch lower(method(1))
  case 'l'
    w = linear(X,v,[xi(:) yi(:) zi(:)]);
  case 'n'
    w = nearest(X,v,[xi(:) yi(:) zi(:)]);
end
w = reshape(w,size(xi));
end


%------------------------------------------------------------
function vi = linear(x,v,xi)
%LINEAR Triangle-based linear interpolation

dt = delaunayTriangulation(x);
if(isreal(v))
    F = scatteredInterpolant(dt.Points,v);
    vi = F(xi);
else
    vre = real(v);
    vim = imag(v);
    F = scatteredInterpolant(dt.Points,vre);
    vire = F(xi);
    F.V = vim;
    viim = F(xi);
    vi = complex(vire,viim);
end
end

%------------------------------------------------------------
function vi = nearest(x,v,xi)
%NEAREST Triangle-based nearest neightbor interpolation

dt = delaunayTriangulation(x);
k = dt.nearestNeighbor(xi);
vi = k;
d = find(isfinite(k));
vi(d) = v(k(d));
end