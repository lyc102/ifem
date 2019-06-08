function [p,t] = circlemesh(x,y,r,h)
%% CIRCLEMESH generates a mesh for a circle
%
% [node,elem] = circlemesh(x,y,r,h) generates a quasi-uniform mesh with
% size h of the circle centered at (x,y) with readius r.
%
% Example
%
%   [node,elem] = circlemesh(0,0,1,0.1);
%   showmesh(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

[p,t] = odtmesh2d(@fd,@huniform,h,[-1,-1;1,1],[],0,x,y,r);

    function s = fd(p,x,y,r)
    s = sqrt(sum((p(:,1)-x).^2+(p(:,2)-y).^2,2))-r;
    end

    function h = huniform(p,varargin)
    h = ones(size(p,1),1);
    end
end