function [node,elem,HB] = cubemesh(box,h)
%% CUBEMESH a uniform mesh of a cube
%
% [node,elem,HB] = cubemesh([x0,x1,y0,y1,z0,z1],h) enerates a uniform mesh
% of the cube [x0,x1]*[y0,y1]*[z0,z1] with mesh size h.
%
% Example
%
%  [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
%  showmesh3(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

x0 = box(1); x1 = box(2); 
y0 = box(3); y1 = box(4);
z0 = box(5); z1 = box(6);
node = [x0,y0,z0; x1,y0,z0; x1,y1,z0; x0,y1,z0; ...
        x0,y0,z1; x1,y0,z1; x1,y1,z1; x0,y1,z1];
elem = [1 2 3 7; 1 4 3 7; 1 5 6 7; 1 5 8 7; 1 2 6 7; 1 4 8 7];
n = ceil(log2(abs(x1-x0)/h));
for k = 1:n
    [node,elem] = uniformrefine3(node,elem);  
end
% Set this as an initial grid
N0 = size(node,1);
HB = zeros(N0,4);
HB(1:N0,1:3) = repmat((1:N0)',1,3); 