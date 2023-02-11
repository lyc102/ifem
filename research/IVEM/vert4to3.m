function [vert3,id,dist] = vert4to3(vert4)

%% USAGE: choose three points among four points such that the plane spanned
%         by the chosen three points has the shortest distance to the 
%         remaining point 

% INPUTS: 
% vert4 --- 4-by 3 matrix stores the coordinates of four points.

% OUTPUTS:
% vert3 --- 3-by 3 matrix stores the coordinates of three chosen points.
% id --- the index of points to be deleted.
% dist --- 4-by-1 vector. dist(i) is the distance from Ai to the plane
%          spanned the the rest.

%% 
dist = zeros(1,4);
for i = 1:4
    vert = vert4; vert(i,:) = [];
    n = cross(vert(2,:)-vert(1,:),vert(3,:)-vert(1,:));
    n = n/norm(n);
    L = @(x,y,z) dot([x-vert(1,1),y-vert(1,2),z-vert(1,3)],n);
    dist(i) = abs(L(vert4(i,1),vert4(i,2),vert4(i,3)));
end
[~,id] = min(dist);
vert3 = vert4;
vert3(id,:) = [];
