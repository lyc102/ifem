function mesh = genMesh3D(domain, nx, ny, nz)

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
c1 = [1,2,3,7]; 
c2 = [1,6,2,7]; 
c3 = [1,5,6,7];
c4 = [1,8,5,7];
c5 = [1,4,8,7];
c6 = [1,3,4,7];
t = [T(:,c1); T(:,c2); T(:,c3); T(:,c4); T(:,c5); T(:,c6)];
mesh = struct('p',p,'t',t,'T',T);