function [node,elem,bdFlag] = uniformrefinehex(node,elem,bdFlag)
%% UNIFORMREFINEQUAD uniformly refine a 3-D hex mesh.
%
% [node,elem] = uniformrefinehex(node,elem) divides each hex into 8 small
% hex by connecting the middle points of the opposited edge of every quad.
% 
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure
T = auxstructurehex(elem);
face = double(T.face);
elem2face = double(T.elem2face);
clear T;

%TODO: finish it
