function [node,elem] = cubehexmesh(cube,h)
%% CUBEHEXMESH uniform mesh of cube
%
% [node,elem] = cubehexmesh([x0,x1,y0,y1,z0,z1],h) generates a uniform
% cubic mesh of the cube [x0,x1]*[y0,y1]*[z0,z1] with mesh size h.
%
% Example
%  [node,elem] = cubehexmesh([0 1 0 2 0 3],0.5);
%
% See also: squarequadmesh, cubehexgradmesh
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

x0 = cube(1); x1 = cube(2); 
y0 = cube(3); y1 = cube(4);
z0 = cube(5); z1 = cube(6);
[z,x,y] = ndgrid(z0:h:z1,x0:h:x1,y0:h:y1);
node = [x(:),y(:),z(:)];
% findnode3(node); view([83 10]);

%% Generate elements
ni = size(x,1); % number of rows
nj = size(x,2); % number of columns
nk = size(x,3); % number of pages
nij = ni*nj;    % number of elements in one page
N = size(node,1);
nodeidx = reshape(1:N,ni,nj,nk);
t2nidxMap = nodeidx(1:ni-1,1:nj-1,1:nk-1);
s = t2nidxMap(:);
elem = [s s+ni s+ni+nij s+nij s+1 s+ni+1 s+ni+nij+1 s+nij+1]; 
% s+1 s+ni*nj
%  | /
%  s - s + ni 
%
%       6 --- 8
%      /|    /|
%     2 --- 4 |
%     | 5 --| 7
%     |/    |/
%     1 --  3