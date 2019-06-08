function [node,elem,elem2dof] = squaremeshp(square,h,p)
%% SQUAREMESHP uniform mesh of a square with periodic condition
%
% [node,elem,elem2dof] = squaremeshp([x0,x1,y0,y1],h,'x') generates a
% uniform mesh of the square [x0,x1]*[y0,y1] with mesh size h with periodic
% condition on x-direction.
%
% [node,elem,elem2dof] = squaremeshp([x0,x1,y0,y1],h,'y') generates a
% uniform mesh of the square [x0,x1]*[y0,y1] with mesh size h with periodic
% condition on y-direction.
%
% [node,elem,elem2dof] = squaremeshp([x0,x1,y0,y1],h,'xy') generates a
% uniform mesh of the square [x0,x1]*[y0,y1] with mesh size h with periodic
% condition on both x and y-direction.
%
% Example
%
%   [node,elem,elem2dof] = squaremeshp([0,1,0,1],0.2,'xy');
%   set(gcf,'Units','normal'); 
%   set(gcf,'Position',[0,0,0.6,0.4]);
%   subplot(1,2,1); showmesh(node,elem); findnode(node);
%   subplot(1,2,2); showmesh(node,elem); findnodedof(node,elem,elem2dof);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

[node,elem] = squaremesh(square,h);
N = size(node,1);  % number of vertices
n = sqrt(N);       % number of vertices in one side
nodeMap = (1:N)';
indexMatrix = reshape(1:n^2,n,n);
switch upper(p)
    case 'X'
        nodeMap(indexMatrix(:,n)) = indexMatrix(:,1);
        elem2dof = nodeMap(elem);
    case 'Y'
        nodeMap(indexMatrix(n,:)) = indexMatrix(1,:);
        indexMap = zeros(n^2,1);
        coarseidx = indexMatrix(1:n-1,1:n);
        indexMap(coarseidx) = 1:(n-1)*n;
        elem2dof = nodeMap(elem);
        elem2dof = indexMap(elem2dof);
    case 'XY'
        nodeMap(indexMatrix(:,n)) = indexMatrix(:,1);
        nodeMap(indexMatrix(n,:)) = indexMatrix(1,:);
        indexMap = zeros(n^2,1);
        coarseidx = indexMatrix(1:n-1,1:n-1);
        indexMap(coarseidx) = 1:(n-1)^2;
        indexMap(indexMatrix(1,n)) = 1;
        elem2dof = nodeMap(elem);
        elem2dof = indexMap(elem2dof);
end

