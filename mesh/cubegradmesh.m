function [node,elem,HB] = cubegradmesh(xh,yh,zh)
%% CUBEGRADMESH a uniform mesh of a cube
%
%
% Example
%     xh = 0:0.125:1;
%     yh = 0:0.125:1;
%     zh = 0:0.125:0.5;
%     [node,elem] = cubegradmesh(xh,yh,zh);
%     showmesh3(node,elem);
%
%   xh = (0:0.125:1).^2;
%   yh = (0:0.125:1).^2;
%   zh = 0:0.125:0.5;
%   [node,elem] = cubegradmesh(xh,yh,zh);
%   showmesh3(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 


%       8 --- 7
%      /|    /|
%     5 --- 6 |      z
%     | 4 --| 3      |  y
%     |/    |/       | /
%     1 --  2        o --- x
%
% The order of vertices is important for the uniform refinement.

[x,y,z] = ndgrid(xh,yh,zh);
node = [x(:),y(:),z(:)];
[nx,ny,nz] = size(x);
elem = zeros(6*(nx-1)*(ny-1)*(nz-1),4);
indexMap = reshape(1:nx*ny*nz,nx,ny,nz);
localIndex = zeros(8,1);
idx = 1;
for k = 1:nz-1
    for j = 1:ny-1
        for i = 1:nx-1
            localIndex(1) = indexMap(i,j,k);
            localIndex(2) = indexMap(i+1,j,k);
            localIndex(3) = indexMap(i+1,j+1,k);
            localIndex(4) = indexMap(i,j+1,k);
            localIndex(5) = indexMap(i,j,k+1);
            localIndex(6) = indexMap(i+1,j,k+1);
            localIndex(7) = indexMap(i+1,j+1,k+1);
            localIndex(8) = indexMap(i,j+1,k+1);
            elem(idx:idx+5,:) = localIndex([1 2 3 7; 1 4 3 7; 1 5 6 7;...
                                            1 5 8 7; 1 2 6 7; 1 4 8 7]);
            idx = idx + 6;
        end
    end
end
% Set this as an initial grid
N0 = size(node,1);
HB = zeros(N0,4);
HB(1:N0,1:3) = repmat((1:N0)',1,3);