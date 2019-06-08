function [node,elem] = cuboidmesh(cuboid,h)
%% CUBEHEXMESH uniform mesh of cuboid
%
% [node,elem] = cuboidmesh([x0,x1,y0,y1,z0,z1],h) generates a uniform mesh of the
% cuboid [x0,x1]*[y0,y1]*[z0,z1] with mesh size h.
%
% Example
%
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

x0 = cuboid(1); x1 = cuboid(2); 
y0 = cuboid(3); y1 = cuboid(4);
z0 = cuboid(5); z1 = cuboid(6);

[x,y,z] = meshgrid(x0:h:x1,y0:h:y1,z0:h:z1);

node = [x(:),y(:),z(:)];

%% Generate elements
nx = size(x,1) - 1; % number of cells in x-direction
ny = size(y,2) - 1; % number of cells in y-direction
nz = size(z,2) - 1; % number of cells in z-direction

elem = zeros(nx*ny*nz,8);
cellidx = 1:nx*ny*nz;

[i, j, k] = ind2sub([nx ny nz],cellidx);  % index of cells in subscript form

s =[ nx+1, ny+1, nz+1];

elem(cellidx,1) = sub2ind3(s,i,j,k);
elem(cellidx,2) = sub2ind3(s,i,j+1,k);
elem(cellidx,3) = sub2ind3(s,i+1,j+1,k);
elem(cellidx,4) = sub2ind3(s,i+1,j,k);
elem(cellidx,5) = sub2ind3(s,i,j,k+1);
elem(cellidx,6) = sub2ind3(s,i,j+1,k+1);
elem(cellidx,7) = sub2ind3(s,i+1,j+1,k+1);
elem(cellidx,8) = sub2ind3(s,i+1,j,k+1);

    function idx = sub2ind3(siz,i,j,k)
        nr = siz(1); nc = siz(2); nv = siz(3);
        if (max(j)>nc) || (max(i)>nr) || (max(k)>nv)
            error(message('MATLAB:mysub2ind:IndexOutOfRange'));
        end
        idx = (k-1)*nr*nc + (j-1)*nr+i;
    end



end
