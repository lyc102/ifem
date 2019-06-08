function [node,elem] = squaremesh(square,h)
%% SQUAREMESH uniform mesh of a square
%
% [node,elem] = squaremesh([x0,x1,y0,y1],h) generates a uniform mesh of the
% rectangle [x0,x1]*[y0,y1] with mesh size h.
%
% Example
%
%   [node,elem] = squaremesh([0,1,0,1],0.2);
%   showmesh(node,elem);
%   findnode(node);
%   findelem(node,elem);
%
% See also: squarequadmesh, cubehexmesh
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Generate nodes
x0 = square(1); x1 = square(2); 
y0 = square(3); y1 = square(4);
[x,y] = meshgrid(x0:h:x1,y0:h:y1);
node = [x(:),y(:)];

%% Generate elements
ni = size(x,1); % number of rows
N = size(node,1);
t2nidxMap = 1:N-ni;
topNode = ni:ni:N-ni;
t2nidxMap(topNode) = [];
k = (t2nidxMap)';
elem = [k+ni k+ni+1 k; k+1 k k+ni+1];
% 4 k+1 --- k+ni+1 3  
%    |        |
% 1  k  ---  k+ni  2