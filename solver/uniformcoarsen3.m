function [elem,HB,newHB,tree]= uniformcoarsen3(elem,HB)
%% UNIFORMCOARSEN3 uniform coarsening in 3-D
%
% [elem,HB,newHB] = uniformcarsen3(elem,HB) remove all good-to-coarsen nodes.
% It is mainly used to get multilevel decomposition in multigrid methods.
% The input matrix elem stands for the fine mesh and the output one for the
% coarse mesh, whose nodal indices are shifted and thus different with the
% fine mesh. 
% 
% Unlike the two dimensional case, an additional matrix HB is introduced.
% The local index of HB is 2 -- 1 -- 3, i.e., HB(:,1) is the middle point
% of the edge formed by HB(:,2:3). HB(:,4) is used to store the generation
% of the node HB(:,1).
% 
% In the output, newHB(:,1) are all removed nodes and newHB(:,2:3) are two
% neighboring nodes whose indices are in the coarse mesh. 
%
% See also: coarsen3, bisect3, uniformcoarsen, MGP1, MGP2, MG3P1, MG3P2
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Find good-to-coarsen nodes
N = max(elem(:));	           
NT = size(elem,1);             
generation = zeros(N,1);      
generation(HB(:,1)) = HB(:,4);     
valence = accumarray(elem(:),ones(4*NT,1),[N 1]);
valenceNew = accumarray(elem(:,4),ones(NT,1), [N 1]); % for newest nodes only
isGoodNode = (valence == valenceNew) & (generation > 0);
if ~any(isGoodNode)
    newHB = []; tree = [];
    return
else
    nGoodNode = sum(isGoodNode);     
    newHB = zeros(nGoodNode,3);      
    newHB(:,1) = find(isGoodNode);   
    newHB(:,2) = HB(newHB(:,1),2);   
    newHB(:,3) = HB(newHB(:,1),3);
end

%% Remove good-to-coarsen nodes
t = find(isGoodNode(elem(:,4)));     % elements containing good nodes
leftNode = HB(elem(t,4),2);          % element-wise left node of good nodes
idx = (elem(t,1)==leftNode) | (elem(t,2)==leftNode) | (elem(t,3)==leftNode);
tl = t(idx);                     
tr = t(~idx);  
Nr = length(tl); % number of refined elements
tree = zeros(Nr,3);
% sort tl and tr by the common vertices such that tl and tr matches
if nargout == 4
    leftNode = HB(elem(tl,4),2);
    rightNode = HB(elem(tr,4),3);
    temptl = zeros(Nr,2);   
    idx = (elem(tl,1) == leftNode);
    temptl(idx,:) = elem(tl(idx),[2 3]);
    idx = (elem(tl,2) == leftNode);
    temptl(idx,:) = elem(tl(idx),[1 3]);
    idx = (elem(tl,3) == leftNode);
    temptl(idx,:) = elem(tl(idx),[1 2]);
    temptl = sort(temptl,2);
    temptr = zeros(Nr,2);
    idx = (elem(tr,1) == rightNode);
    temptr(idx,:) = elem(tr(idx),[2 3]);
    idx = (elem(tr,2) == rightNode);
    temptr(idx,:) = elem(tr(idx),[1 3]);
    idx = (elem(tr,3) == rightNode);
    temptr(idx,:) = elem(tr(idx),[1 2]);
    temptr = sort(temptr,2);
    [temptl, Il] = sortrows([temptl elem(tl,4)]);
    [temptr, Ir] = sortrows([temptr elem(tr,4)]);
    tl = tl(Il);
    tr = tr(Ir);
    tree(:,3) = tr;
    tree(:,1) = tl;
    tree(:,2) = tl;
end
% coarsen tl to t
elem(tl,4) = HB(elem(tl,4),3);       % replace new node by right node

%% Sort element nodes by generations
% the newest node is the node with maxmum generation
if (length(tl)==1)    
    idx = max(transpose(generation(elem(t1,:)))); 
else
    [tempvar,idx] = max(generation(elem(tl,:)),[],2);  %#ok<*ASGLU>
end
elem(tl((idx==1)),1:4) = elem(tl((idx==1)),[2 4 3 1]);
elem(tl((idx==2)),1:4) = elem(tl((idx==2)),[3 4 1 2]);
elem(tl((idx==3)),1:4) = elem(tl((idx==3)),[4 2 1 3]);

%% Clean and shift index
elem(tr,:) = [];                      
inCoarse = true(NT,1);
inCoarse(tr) = false;
elemidxMap = zeros(NT,1);
elemidxMap(inCoarse) = 1:size(elem,1);
if nargout == 4
    tree(:,1) = elemidxMap(tree(:,1));
end
HB(isGoodNode,:) = [];                
indexMap = zeros(N,1);                
indexMap(~isGoodNode)= 1:N-nGoodNode; 
elem = indexMap(elem);                
HB(:,1:3) = indexMap(HB(:,1:3));      
newHB(:,2:3) = indexMap(newHB(:,2:3));