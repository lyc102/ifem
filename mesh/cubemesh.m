function [node,elem,HB] = cubemesh(box,varargin)
%% CUBEMESH a uniform mesh of a cube
%
% [node,elem,HB] = cubemesh([x0,x1,y0,y1,z0,z1],h) generates a uniform mesh
% of the cube [x0,x1]*[y0,y1]*[z0,z1] with certain mesh size h. 
%
% Default mode: the number of cuts along x, y, z direction are equal.
% Additional option: different mesh size along x, y, z direction.
%
% Example 1
%
%  [node,elem] = cubemesh([-1,1,-1,1,-1,1]);
%  showmesh3(node,elem);
%
% Example 2
%
%  [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
%  showmesh3(node,elem);
% 
% Example 3
%
%  option.hx = 1/2; option.hy = 1/3; option.hz = 1;
%  [node,elem,HB] = cubemesh([0,1,0,1,0,1],option);
%  showmesh3(node,elem);
%
% The ordering of vertices is important for the uniform refinement. See
% uniformrefine3.
% 
% See also:
%    uniformrefine3
%
% Different mesh size in x-y-z directions are added by Shuhao Cao.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

x0 = box(1); x1 = box(2); 
y0 = box(3); y1 = box(4);
z0 = box(5); z1 = box(6);

node = [x0,y0,z0; x1,y0,z0; x1,y1,z0; x0,y1,z0; ...
        x0,y0,z1; x1,y0,z1; x1,y1,z1; x0,y1,z1];
elem = [1 2 3 7; 1 4 3 7; 1 5 6 7; 1 5 8 7; 1 2 6 7; 1 4 8 7];
% The order of vertices is important for the uniform refinement.

if nargin ==1
    h = min(abs([x1- x0, y1-y0, z1-z0]/2));
    
    n = max(ceil(log2(abs([x1- x0, y1-y0, z1-z0]/h))));
    for k = 1:n
        [node,elem] = uniformrefine3(node,elem);
    end
    
elseif (nargin >= 2) && isnumeric(varargin{1})
    h = varargin{1};
    n = max(ceil(log2(abs([x1- x0, y1-y0, z1-z0]/h))));
    for k = 1:n
        [node,elem] = uniformrefine3(node,elem);
    end
      
elseif  (nargin >= 2) && ~isnumeric(varargin{1})
    option = varargin{1};
    hx = option.hx; hy = option.hy; hz = option.hz;
    
    hx = (x1- x0)/floor(abs((x1- x0)/hx));
    hy = (x1- x0)/floor(abs((x1- x0)/hy));
    hz = (x1- x0)/floor(abs((x1- x0)/hz));
  
    [x,y,z] = ndgrid(x0:hx:x1,y0:hy:y1,z0:hz:z1);
    node = [x(:),y(:),z(:)];

    matlabversion = version();
    if str2double(matlabversion(end-5:end-2)) <= 2013
        T = DelaunayTri(x(:),y(:),z(:)); %#ok<*DDELTRI>
        elem = T.Triangulation;
    else
        T = delaunayTriangulation(x(:),y(:),z(:));
        elem = T.ConnectivityList;
    end
    elem = label3(node,elem);
end

% Set this as an initial grid
N0 = size(node,1);
HB = zeros(N0,4);
HB(1:N0,1:3) = repmat((1:N0)',1,3);