function [node,elem] = optsurfacemesh(node,elem,surfacedata, gamma, maxit)
%% OPTSURFACEMESH optimise surface mesh by lloyd method.
%
%
% Author: Huayi Wei <huayiwei1984@gmail.com>
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if nargin <= 4
    maxit = 20;
end

if nargin <= 3
    gamma = 1;
end

if nargin < 3
   fprintf('I need three args: node,elem and surfacedata!');
   return
end

for i = 1: maxit
    node = masscenter(node,elem,surfacedata,gamma);
    node = surfacedata.project(node);
    flag = 1;
    while flag
        [elem,flag] = edgeswap(node,elem);
    end
end
