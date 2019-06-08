function [Pro,freeEdgec] = transferedge(elemc,elemf,tree,freeEdge)
%% TRANSFEREDGE transfer operator between bisected grids
%
%

NTc = size(elemc,1);

if ~exist('freeEdge','var'),freeEdge = []; end

%% First coarsen
isMarkedElem([tree(:,2); tree(:,3)]) = true;
[tempnode,tempelem,~,~,temptree] = coarsen([],elemf,isMarkedElem);
[Pro_1,tempFreeEdge] = transferedgecoarsen(tempelem,elemf,temptree,freeEdge);

%% Index map between fine mesh and temp mesh
NTf = size(elemf,1); NTt = size(tempelem,1);
inTemp = true(NTf,1);
inTemp(temptree(:,3)) = false;  
elemf2t = zeros(NTf,1);    
elemf2t(inTemp) = 1:NTt; % fine to temp index map of elements 

%% Second coarsen
idx = (tree(:,1) <= NTc);
isMarkedElem([tree(idx,2); tree(idx,3)]) = true; % idx in elemf
markedElem = elemf2t(isMarkedElem); % idx in tempelem
markedElem(markedElem == 0) = [];
[tempnode,elemc,~,~,treec] = coarsen([],tempelem,markedElem);
[Pro_2,freeEdgec] = transferedgecoarsen(elemc,tempelem,treec,tempFreeEdge);

%% Prolongation
Pro = Pro_1*Pro_2;