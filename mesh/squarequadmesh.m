function [node,elem,T] = squarequadmesh(square,h)
%% SQUAREQUADMESH uniform rectangular mesh of a square
%
% [node,elem] = squarequadmesh([x0,x1,y0,y1],h) generates a uniform mesh of the
% rectangle [x0,x1]*[y0,y1] with mesh size h.
%
% [node,elem,T] = squarequadmesh([x0,x1,y0,y1],h) with additional output T
% contains auxiliary structure: edge, elem2edge, edge2elem, bdEdge
%
% Example
%
%   [node,elem] = squarequadmesh([0,1,0,2],0.5);
%   showmesh(node,elem);
%   findnode(node);
%   findquadelem(node,elem);
%
% See also: squaremesh, squaregradmeshquad, cubehexmesh
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Generate nodes
x0 = square(1); x1 = square(2); 
y0 = square(3); y1 = square(4);
% [x,y] = meshgrid(x0:h:x1,y1:-h:y0);
[x,y] = ndgrid(x0:h:x1,y0:h:y1);
node = [x(:),y(:)];

%% Generate elements
ni = size(x,1); % number of rows
nj = size(x,2);
N = size(node,1);
nodeidx = reshape(1:N,ni,nj);
t2nidxMap = nodeidx(1:ni-1,1:nj-1);
k = t2nidxMap(:);
elem = [k k+ni k+ni+1 k+1];
NT = size(elem,1);
% 4 k+1 --- k+ni+1 3  
%    |        |
% 1  k  ---  k+ni  2

%% Check if further structure is needed
if nargout<=2
    return;
end

%% Generate edges
Ne1 = (ni-1)*nj;
Ne2 = ni*(nj-1);
edge = zeros(Ne1+Ne2,2,'int32');
% vertical edges
e2nidxMap = nodeidx(1:ni-1,1:nj);
k = e2nidxMap(:);
edge(1:Ne1,:) = [k k+1];
% horizontal edges
e2nidxMap = nodeidx(1:ni,1:nj-1);
k = e2nidxMap(:);
edge(Ne1+(1:Ne2),:) = [k k+ni];

%% Generate index map between element and edges
elem2edge = zeros(NT,4,'int32');
edge2elem = zeros(Ne1+Ne2,2,'int32');
k = 1:NT;
veidx = k; % veidx = sub2ind([ni-1,nj],i,j);
heidx = Ne1 + (1:Ne2); % linear indexing
heidx = reshape(heidx,ni,nj-1); % matrix indexing
heidx = heidx(1:ni-1,:); % without top edges
elem2edge(k,1) = heidx(:);
elem2edge(k,3) = heidx(:)+1;
elem2edge(k,2) = veidx + (ni-1);
elem2edge(k,4) = veidx;
% vertical edges
edge2elem(veidx,2) = k;   % right element
edge2elem(veidx+ni-1,1) = k; % left element
% horizontal edges
edge2elem(heidx,1) = k;
edge2elem(heidx+1,2) = k;

%% Generate boundary edges
% for boundary edges, set the another element idx to zero
idx1 = edge2elem(:,1) == 0; % left and top
idx2 = edge2elem(:,2) == 0;
bdEdge = [edge(idx1,[2 1]); edge(idx2,:)];  
% switch edge index on the left and top to have a consistent orientation

%% Auxstructure
T = struct('edge',edge,'elem2edge',elem2edge,'edge2elem',edge2elem,'bdEdge',bdEdge);
% neighbor and bdElem can be easily get from 2-D index of all elements