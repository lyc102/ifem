%% FLOWPASTCYLINDER Generate a rectangular mesh with a circular hole
%
% It is used in the benchmark problem in CFD: flow past a cylinder. The
% domain is (0,2.2)*(0,0.41) - {(x,y) | (x-0.2)^2+(y-0.2)^2<=0.05^2}. The
% mesh is refined at the circle boundary.
%
% Example
%   flowpastcylinder;
%   showmesh(node,elem);
%
% See also squaremesh, circlemesh
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

fd = inline('ddiff(drectangle(p,0,2.2,0,0.41),dcircle(p,0.2,0.2,0.05))','p');
fh = inline('min(sqrt(sum((p-0.2).^2,2))-(0.05-0.0125),0.346)','p');
[node,elem] = odtmesh2d(fd,fh,0.0035,[0,0;2.2,0.41],[0,0;2.2,0;2.2,0.41;0,0.41]);