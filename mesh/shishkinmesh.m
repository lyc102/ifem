function [node,elem] = shishkinmesh(square,side,epsilon,n,c)
%% SHISHKINMESH layer-adapted piecewise uniform mesh
%
% [node,elem] = shishkinmesh(square,'L',epsilon,n,c) generates a piecewise
% uniform mesh of the square [x0,x1]*[y0,y1] refined near the left boundary
% of the square. epsilon is the length of the boundary layer and n is the
% number vertices in one direction. The parameter c is used in the
% computation of the transition point tau = c*epsilon*log(n)*(x1-x0). It is
% optional and the default value is 2. 
%
% The second input argument is a string to specify the location of the
% boundary layer and now can be 'L','R','U','B'.
%
% Example
%   [node,elem] = shishkinmesh([0 1 0 1],'U',0.01,10);
%   subplot(1,2,1); showmesh(node,elem);
%   [node,elem] = shishkinmesh([0 1 0 1],'B',0.01,10);
%   subplot(1,2,2); showmesh(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

x0 = square(1); x1 = square(2); 
y0 = square(3); y1 = square(4);
h = (x1-x0)/n;     
[node,elem] = squaremesh(square,h);
N = size(node,1);  % number of vertices
n = sqrt(N);       % number of vertices in one side
midn = ceil(n/2);
if ~exist('c','var'), c = 2; end
switch upper(side)
    case 'L'
        tau = c*epsilon*log(n)*(x1-x0);
        if tau>1/2, return; end
        h = tau/(midn-1);
        H = (x1-x0-tau)/(n - midn);
        scalex = repmat(x0:h:x0+tau,n,1);
        node(1:midn*n,1) = scalex(:);
        scalex = repmat(x0+tau+H:H:x1,n,1);
        node(midn*n+1:n^2,1) = scalex(:);
    case 'R'
        tau = c*epsilon*log(n)*(x1-x0);
        if tau>1/2, return; end
        h = tau/(n - midn);
        H = (x1-x0-tau)/(midn-1);
        scalex = repmat(x0:H:(x1-tau),n,1);
        node(1:midn*n,1) = scalex(:);
        scalex = repmat((x1-tau)+h:h:x1,n,1);
        node(midn*n+1:n^2,1) = scalex(:);
    case 'U'
        tau = c*epsilon*log(n)*(y1-y0);
        if tau>1/2, return; end
        h = tau/(n - midn);
        H = (y1-y0-tau)/(midn-1);
        scaley = repmat((y0:H:(y1-tau))',1,n);
        idx = repmat((0:midn-1)',1,n) + repmat((0:n-1)*n+1,midn,1);
        node(idx,2) = scaley(:);
        scaley = repmat(((y1-tau)+h:h:y1)',1,n);
        idx = repmat((midn:n-1)',1,n) + repmat((0:n-1)*n+1,n-midn,1);
        node(idx,2) = scaley(:);
    case 'B'
        tau = c*epsilon*log(n)*(y1-y0);
        h = tau/(midn-1);
        H = (y1-y0-tau)/(n - midn);
        scaley = repmat((y0:h:y0+tau)',1,n);
        idx = repmat((0:midn-1)',1,n) + repmat((0:n-1)*n+1,midn,1);
        node(idx(:),2) = scaley(:);
        scaley = repmat((y0+tau+H:H:y1)',1,n);
        idx = repmat((midn:n-1)',1,n) + repmat((0:n-1)*n+1,n-midn,1);
        node(idx,2) = scaley(:);
end