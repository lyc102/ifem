function [node,elem,Ty] = squaregradmeshquad(square,h,gamma)
%% SQUAREGRADMESHQUAD 
%  Generate a graded mesh in y direction use tranformation y^gamma.
%
%   [node,elem] = squaregradmesh([0,1,0,2],0.1,0.2);
%   showmesh(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Change coordinate of y
x0 = square(1); x1 = square(2); 
y0 = square(3); y1 = square(4);
Tx = x0:h:x1; % uniform in x
My = (y1-y0)/h;
k = 0:My;
% Ty = (k/My).^gamma*(y1-y0) + y0;
Ty = gradmap(y0,y1,gamma,h);
[x,y] = meshgrid(Tx,Ty);
node = [x(:),y(:)];

%% Generate elements
ni = size(x,1); % number of rows
nj = size(x,2);
N = size(node,1);
nodeidx = reshape(1:N,ni,nj);
t2nidxMap = nodeidx(1:ni-1,1:nj-1);
k = t2nidxMap(:);
elem = [k k+ni k+ni+1 k+1];