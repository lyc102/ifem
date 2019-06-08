function cubeIdx = getCubeIdx(p, box, h)
%% GETCUBEIDX get the index of small cube containing p
%
%
%
%

ni = floor((box(2) - box(1))/h)+1;
nj = floor((box(4) - box(3))/h)+1;

x = p(:, 1) - box(1);
y = p(:, 2) - box(3);
z = p(:, 3) - box(5);

i = floor(x/h)+1;
j = floor(y/h)+1;
k = floor(z/h)+1;

nij = (ni-1)*(nj-1);

cubeIdx = (j-1)*nij + (i-1)*(ni-1) + k;