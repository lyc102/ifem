function [node,elem] = cvtuniformmesh(n)
% n is the number of intervales in one coordinate direction

h = 2/n;
[x y z] = meshgrid(-1:h:1,-1:h:1,-1:h:1);
[cx cy cz] = meshgrid(-1+h/2:h:1-h/2,-1+h/2:h:1-h/2,-1+h/2:h:1-h/2);
node(:,1) = [x(:); cx(:)];
node(:,2) = [y(:); cy(:)];
node(:,3) = [z(:); cz(:)];
elem = delaunayn(node);
elem = fixorder3(node,elem);