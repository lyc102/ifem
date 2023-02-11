function mesh = genMesh3DType2(domain, nx, ny, nz)

%% Usage: mesh structure of a uniform cubic partition
%
% INPUTS:
% domain --- cubic domain = [xmin, xmax, ymin, ymax, zmin, zmax].
% nx --- the number of uniform partition in x direction.
% ny --- the number of uniform partition in y direction.
% nz --- the number of uniform partition in z direction.
%
% option.meshinfo --- basic (default): generate only p and t.
%                     all: enriched mesh information
%
% OUTPUTS:
% mesh --- a struct data contains mesh information.
% 
% Last Modified: 08/07/2020 by Xu Zhang
%
%                A8-------------------A7        The Cube is divided into
%                /|                   /|        six congruent tetrahedrons
%               / |                  / |        
%              /  |                 /  |        
%             /   |                /   |        (1) A1-A2-A3-A7
%            /    |               /    |        (2) A1-A6-A2-A7
%          A5-----+-------------A6     |        (3) A1-A5-A6-A7
%           |     |             |      |        (4) A1-A8-A5-A7
%           |     |             |      |        (5) A1-A4-A8-A7
%           |     |             |      |        (6) A1-A3-A4-A7
%           |     A4------------+------A3
%           |     /             |      /
%           |    /              |     /
%           |   /               |    /
%           |  /                |   /
%           | /                 |  /
%           |/                  | /
%          A1-------------------A2
%
%% 1. Generate basic mesh info: p t
[p,T] = genMesh3DRectPT(domain, nx, ny, nz);
c1 = zeros(5,4);
c1(1,:) = [1,2,4,5]; 
c1(2,:) = [2,7,5,6]; 
c1(3,:) = [3,2,4,7];
c1(4,:) = [5,7,4,8];
c1(5,:) = [2,4,5,7];
c2 = zeros(5,4);
c2(1,:) = [1,6,8,5]; 
c2(2,:) = [1,3,6,2]; 
c2(3,:) = [8,6,3,7];
c2(4,:) = [1,8,3,4];
c2(5,:) = [1,3,8,6];

CBxy = zeros(nx,ny);
CBxy(1:2:nx,2:2:ny)=1;
CBxy(2:2:nx,1:2:ny)=1;
CB = repmat(CBxy,1,1,nz);
CB(:,:,1:2:nz) = repmat(1-CBxy,1,1,length(1:2:nz));
CB = reshape(CB,[],1);
tid1 = find(CB==1); tid2 = find(CB==0);

t = T(:,c1(1,:));
t(tid2,:) = T(tid2,c2(1,:));
for i = 2:5
    tmp = T(:,c1(i,:));
    tmp(tid2,:) = T(tid2,c2(i,:));
    t = [t;tmp];
end
  
mesh = struct('p',p,'t',t,'T',T);