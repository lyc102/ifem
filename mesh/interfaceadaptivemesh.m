function [node,elem,phiValue,oldNodeIdx,oldNode,oldElemIdx,oldElem] = interfaceadaptivemesh(node,elem,phi)
% INTERFACEMESH bisect the initial mesh to get a mesh fiting the interface
% 
%   [node,elem] = INTERFACEADAPTIVEMESH(node,elem,phi,a), where 
%    - node and elem represent an initial mesh enclosing the interface; 
%    - phi is a level set function defining the interface as follows
%       {p | phi(p)=0} represents the interface
%       {p | phi(p) < 0} represents the interior domain enclosed by the interface
%       {p | phi(p) > 0} rpepresents the exterior domain of the interface. 
%  
%  Definition
%   Interface element: the element whose vertices are not in the same side
%   of the interface.
%
%   Interface node: the vertices of all the interface elements. Split the
%   set of interface nodes into two parts: interior interface nodes (phi(p)<=0) and
%   exterior interface nodes (phi(p)>=0);
%
%  Example
%   clear all; close all;
%   node = [-1 -1; 1 -1; 1 1; -1 1];
%   elem = [2 3 1; 4 1 3];
%   [node,elem] = interfacemesh(node,elem,@ellipse);
%
%  See also bisect, circle, ellipse.
%
%  Reference: H. Wei, L. Chen, Y. Huang and B. Zheng. Adaptive Mesh
%  Refinement and Superconvergence for Two Dimensional Interface Problems.
%  SIAM Journal on Scientific Computing, 2014.
%
%  Add by Huayi Wei. Based on discussion with Long Chen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

oldNodeIdx = [];
oldNode = [];
oldElemIdx = [];
oldElem = [];

%% Step 1: Refine mesh such that the set of interface elements is non-empty
N0 = size(node,1);
phiValue = phi(node);

if all(phiValue<0)
    error('The initial mesh is in the interior of the domain the interface surrounds, please use a mesh which can cover the interface!')
end

% uniform bisect the mesh, to make sure there exist at least one mesh node
% which is inside of the interface, namely, phi(p)<0.
while all(phiValue>0)
    [node,elem] = uniformbisect(node,elem);
    N1 = size(node,1);
    phiValue(N0+1:N1) = phi(node(N0+1:N1,:)); % append the phi value
    N0 = N1;
end

display('Finish Step 1: uniform bisect all elements');

%% Step 2: estimate the discrete curvature, then bisect elements with big curvature
% find interface elements, which are elements crossed by the interface and
% restrict the computation to the interface elements
[interfaceElemIdx, ~, isInterfaceElem] = markinterface(elem,phiValue);
interfaceElem = elem(interfaceElemIdx,:);
isInterfaceNode = false(size(node,1),1);
isInterfaceNode(interfaceElem(:))=true;
interfaceNT = size(interfaceElem,1);
T = auxstructure(interfaceElem);
edge = double(T.edge);
edge2elem = T.edge2elem;
clear T;

% case: Fig 2(a) in p6 
tmpAngle = ones(interfaceNT,1)*[90,45,45];
angle = accumarray(interfaceElem(:),tmpAngle(:),[size(node,1) 1]);
is360 = (angle == 360);

% case: Fig 2(b) in p6
isInterfaceBEdge = (edge2elem(:,1) == edge2elem(:,2));
interfaceBEdge = edge(isInterfaceBEdge,:);
v2vI = sparse(interfaceBEdge,interfaceBEdge(:,[2,1]),1,size(node,1),size(node,1));
isLinkNode = (v2vI*isInterfaceNode >= 3);

% case: Fig 2(c) in p6
isInteriorNode = phiValue < 0;
turnAngle = zeros(size(node,1),1);
isIn = isInterfaceNode & isInteriorNode & ~is360 & ~isLinkNode;
turnAngle(isIn) = 180 - (360 - angle(isIn));
isOut = isInterfaceNode & ~isInteriorNode & ~is360 & ~isLinkNode;
turnAngle(isOut ) = (360 - angle(isOut)) -180;
v2v = sparse(edge,edge(:,[2,1]),1,size(node,1),size(node,1));
turnAngleTotal = v2v*turnAngle+turnAngle;

% case: size of interface elements should be small enough
ve = node(elem(:,3),:) - node(elem(:,2),:);
L = sqrt(sum(ve.^2,2));
mL = max(L);
isBigSizeElem = ((L > mL*mL) & isInterfaceElem);
% mark nodes for refinement
isBigCurvNode = (is360 | isLinkNode |(abs(turnAngleTotal)> 170));

if any(isBigCurvNode) || any(isBigSizeElem)
    
    bigCurvElem = interfaceElemIdx(max(isBigCurvNode(interfaceElem),[],2));
    bigSizeElem = find(isBigSizeElem);   

    count = 0;
    while  count < 100
        % bisect marked elements
        [node,elem] = bisect(node,elem,[bigCurvElem; bigSizeElem]);
        
        N1 = size(node,1);
        phiValue(N0+1:N1) = phi(node(N0+1:N1,:));
        N0 = N1;
        
        % find interface elements and nodes
        interfaceElemIdx = markinterface(elem,phiValue);
        interfaceElem = elem(interfaceElemIdx,:);
        isInterfaceNode = false(size(node,1),1);
        isInterfaceNode(interfaceElem(:))=true;
        T = auxstructure(interfaceElem);
        edge = double(T.edge);
        edge2elem = T.edge2elem;
        clear T;
                
        % find triangles containing curves with high curvature
        % case: Fig 2(a) in p6 
        interfaceNT = size(interfaceElem,1);
        tmpAngle = ones(interfaceNT,1)*[90,45,45];        
        angle = accumarray(interfaceElem(:),tmpAngle(:),[size(node,1) 1]);
        is360 = angle == 360;
        
        % case: Fig 2(b) in p6
        isInterfaceBEdge = edge2elem(:,1) == edge2elem(:,2);
        interfaceBEdge = edge(isInterfaceBEdge,:);
        v2vI = sparse(interfaceBEdge,interfaceBEdge(:,[2,1]),1,size(node,1),size(node,1));
        isLinkNode = (v2vI*isInterfaceNode >= 3);

        % case: Fig 2(c) in p6
        isInteriorNode = (phiValue < 0);
        turnAngle = zeros(size(node,1),1);
        isIn = isInterfaceNode & isInteriorNode & ~is360 & ~isLinkNode;
        turnAngle(isIn) = 180 - (360 - angle(isIn));
        isOut = isInterfaceNode & ~isInteriorNode & ~is360 & ~isLinkNode;
        turnAngle(isOut ) = (360 - angle(isOut)) -180;
        v2v = sparse(edge,edge(:,[2,1]),1,size(node,1),size(node,1));
        turnAngleTotal = v2v*turnAngle+turnAngle;
        
        % control size of interface elements
        ve = node(interfaceElem(:,3),:) - node(interfaceElem(:,2),:);
        L = sqrt(sum(ve.^2,2));
        bigSizeElem = interfaceElemIdx(L > mL^2);
        
        isBigCurvNode = is360 | isLinkNode |(abs(turnAngleTotal)> 170);
        if ~any(isBigCurvNode) && ~any(isBigSizeElem)
            break;
        end
        
        bigCurvElem = interfaceElemIdx(max(isBigCurvNode(interfaceElem),[],2));
    end
end

display('Finish Step 2: bisect interface with big curvature!');

%% Step 3: Update the type of interface triangles
% Two types of triangles in the bisection grids
%  Type A: two edges parallel to coordinates; see Fig 2.6
%  Type B: the longest edge parallel to x or y coordinate; see Fig 2.7

ve = node(elem(:,3),:) - node(elem(:,2),:);
elemType = (abs(ve(:,1)) > eps & abs(ve(:,2)) > eps); 
[interfaceElemIdx, ~,isInterfaceElem] = markinterface(elem,phiValue);
T = auxstructure(elem);
neighbor = T.neighbor;
clear T;
isMark = false(size(elem,1),1);
isMark(interfaceElemIdx) = true;
isMark(neighbor(interfaceElemIdx,1)) = true;

t = cumsum(ones(size(elem,1),1));
tt = neighbor(:,1);
isTransitionTypeBElem = isInterfaceElem & isInterfaceElem(tt) & neighbor(tt,1) ~= t & ~elemType;

%% Deal with the transition Type B elements, see Fig 2.7 (c) and (d) and conditions (C1) and (C2)
% 
%  case 1: a Type B interface element whose second neighbor is a interface
%          and  Type A element.
idx2 = find(isTransitionTypeBElem & isInterfaceElem(neighbor(:,2)) ...
            & elemType(neighbor(:,2)));
if ~isempty(idx2)
    p1 = node(elem(idx2,1),:);
    p2 = node(elem(idx2,2),:);
    p3 = node(elem(idx2,3),:);
    p4 = node(elem(neighbor(idx2,2),1),:);
    m14 = (p1 + p4)/2;
    m34 = (p3 + p4)/2;
    m23 = (p2 + p3)/2;
    m1 = (m23 + p3)/2;
    m2 = (m23 + p2)/2;
    isNotTB = ((phiValue(elem(idx2,3)) .* phi(m1)) > 0 & ...
               (phi(m14) .* phi(m34) < 0)) | phi(m1).*phi(m2) < 0;
    isTransitionTypeBElem(idx2(isNotTB)) = false;
end

%  case 2: a Type B interface element whose third neighbor is a interface
%          and  Type A element.
idx3 = find(isTransitionTypeBElem & isInterfaceElem(neighbor(:,3)) ...
            & elemType(neighbor(:,3)));
if ~isempty(idx3)
    p1 = node(elem(idx3,1),:);
    p2 = node(elem(idx3,2),:);
    p3 = node(elem(idx3,3),:);
    p4 = node(elem(neighbor(idx3,3),1),:);
    m14 = (p1 + p4)/2;
    m24 = (p2 + p4)/2;
    m23 = (p2 + p3)/2;
    m1 = (m23 + p3)/2;
    m2 = (m23 + p2)/2;
    isNotTB = ((phiValue(elem(idx3,2)).*phi(m2)) > 0 & ...
               (phi(m14).*phi(m24) < 0)) | phi(m1).*phi(m2) < 0;
    isTransitionTypeBElem(idx3(isNotTB)) = false;
end

%  case 3: a Type B interface element whose third neighbor is a interface
%          and  Type B element.
idx4 = find(isTransitionTypeBElem & isInterfaceElem(neighbor(:,2)) ...
          & ~elemType(neighbor(:,2)));
if ~isempty(idx4)
    p2 = node(elem(idx4,2),:);
    p3 = node(elem(idx4,3),:);
    m23 = (p2 + p3)/2;
    m1 = (m23 + p3)/2;
    m2 = (m23 + p2)/2;
    isNotTB =  phi(m1).*phi(m2) < 0;
    isTransitionTypeBElem(idx4(isNotTB)) = false;
end

%  case 4: a Type B interface element whose third neighbor is a interface
%          and  Type B element.
idx5 = find(isTransitionTypeBElem & isInterfaceElem(neighbor(:,3)) ...
          & ~elemType(neighbor(:,3)));
if ~isempty(idx5)
    p2 = node(elem(idx5,2),:);
    p3 = node(elem(idx5,3),:);
    m23 = (p2 + p3)/2;
    m1 = (m23 + p3)/2;
    m2 = (m23 + p2)/2;
    isNotTB =  phi(m1).*phi(m2) < 0;
    isTransitionTypeBElem(idx5(isNotTB)) = false;
end
typeB = (isMark & ~elemType & ~isTransitionTypeBElem);

% bisect some transition type B elements to type A
while any(typeB)
    
    [node,elem,~, ~,tree] = bisect(node,elem,typeB);
    N1 = size(node,1);
    phiValue(N0+1:N1) = phi(node(N0+1:N1,:));
    N0 = N1;
    
    elemType = updatetype(elemType,tree(:,[2,3]));
    
    [interfaceElemIdx, ~,isInterfaceElem] = markinterface(elem,phiValue);
    
    T = auxstructure(elem);
    neighbor = T.neighbor;
    clear T;
    isMark = false(size(elem,1),1);
    isMark(interfaceElemIdx) = true;
    isMark(neighbor(interfaceElemIdx,1)) = true;
    
    t = cumsum(ones(size(elem,1),1));
    tt = neighbor(:,1);
    isTransitionTypeBElem = isInterfaceElem & isInterfaceElem(tt) & ...
                            neighbor(tt,1) ~= t & ~elemType;
    
    %  case 1: a Type B interface element whose second neighbor is a interface
    %          and  Type A element.
    idx2 = find(isTransitionTypeBElem & isInterfaceElem(neighbor(:,2)) & ...
                elemType(neighbor(:,2)));
    if ~isempty(idx2)
        p1 = node(elem(idx2,1),:);
        p2 = node(elem(idx2,2),:);
        p3 = node(elem(idx2,3),:);
        p4 = node(elem(neighbor(idx2,2),1),:);
        m14 = (p1 + p4)/2;
        m34 = (p3 + p4)/2;
        m23 = (p2 + p3)/2;
        m1 = (m23 + p3)/2;
        m2 = (m23 + p2)/2;
        isNotTB = ((phiValue(elem(idx2,3)) .* phi(m1)) > 0 & ...
                   (phi(m14) .* phi(m34) < 0)) | phi(m1).*phi(m2) < 0;
        isTransitionTypeBElem(idx2(isNotTB)) = false;
    end
    %  case 2: a Type B interface element whose third neighbor is a interface
    %          and  Type A element.
    idx3 = find(isTransitionTypeBElem & isInterfaceElem(neighbor(:,3)) & ...
                elemType(neighbor(:,3)));
    if ~isempty(idx3)
        p1 = node(elem(idx3,1),:);
        p2 = node(elem(idx3,2),:);
        p3 = node(elem(idx3,3),:);
        p4 = node(elem(neighbor(idx3,3),1),:);
        m14 = (p1 + p4)/2;
        m24 = (p2 + p4)/2;
        m23 = (p2 + p3)/2;
        m1 = (m23 + p3)/2;
        m2 = (m23 + p2)/2;
        isNotTB = ((phiValue(elem(idx3,2)).*phi(m2)) > 0 & ...
                   (phi(m14).*phi(m24) < 0)) | phi(m1).*phi(m2) < 0;
        isTransitionTypeBElem(idx3(isNotTB)) = false;
    end
    
    %  case 3: a Type B interface element whose second neighbor is a interface
    %          and  Type B element.
    idx4 = find(isTransitionTypeBElem & isInterfaceElem(neighbor(:,2)) & ...
               ~elemType(neighbor(:,2)));
    if ~isempty(idx4)        
        p2 = node(elem(idx4,2),:);
        p3 = node(elem(idx4,3),:);
        m23 = (p2 + p3)/2;
        m1 = (m23 + p3)/2;
        m2 = (m23 + p2)/2;
        isNotTB =  phi(m1).*phi(m2) < 0;
        isTransitionTypeBElem(idx4(isNotTB)) = false;
    end
    %  case 4: a Type B interface element whose third neighbor is a interface
    %          and  Type B element.
    idx5 = find(isTransitionTypeBElem & isInterfaceElem(neighbor(:,3)) & ...
                ~elemType(neighbor(:,3)));
    if ~isempty(idx5)
        p2 = node(elem(idx5,2),:);
        p3 = node(elem(idx5,3),:);
        m23 = (p2 + p3)/2;
        m1 = (m23 + p3)/2;
        m2 = (m23 + p2)/2;
        isNotTB =  phi(m1).*phi(m2) < 0;
        isTransitionTypeBElem(idx5(isNotTB)) = false;
    end
    typeB = (isMark & ~elemType & ~isTransitionTypeBElem);    
end

% find interface elements
interfaceElemIdx= markinterface(elem,phiValue);
interfaceElem = elem(interfaceElemIdx,:);
% find interface nodes, which are the vertices of interface elements.
isInterfaceNode = false(size(node,1),1);
isInterfaceNode(interfaceElem(:))=true;
T = auxstructure(elem);
edge = T.edge;
edge2elem = T.edge2elem;
clear T;
% case: Fig 2.7 (b)
isEdge = edge2elem(:,1) ~= edge2elem(:,2);
isEdge = isEdge & (isInterfaceNode(edge(:,1)) | isInterfaceNode(edge(:,2)) |...
    isInterfaceNode(elem(edge2elem(:,1),1)) | isInterfaceNode(elem(edge2elem(:,2),1)));
isEdge = isEdge &  edge2elem(:,3) == 1 & edge2elem(:,4)==1;
isEdge = isEdge & ~elemType(edge2elem(:,1));
isMark = false(size(elem,1),1);
isMark(edge2elem(isEdge,1)) = true;
isMark(edge2elem(isEdge,2)) = true;
[node,elem,~, ~,tree] = bisect(node,elem,isMark);
N1 = size(node,1);
phiValue(N0+1:N1) = phi(node(N0+1:N1,:));
% N0 = N1;
elemType = updatetype(elemType,tree(:,[2,3]));

display('Finish Step 3: transform the type B of  element to type A!');


%% Step 4: Move some interface nodes onto the interface
% 
%  Here we just move interface nodes along the shortest edges of all the
%  type A interface elements. For a such kind of edge, we get the
%  intersection points and compute the distances between the two endpoints
%  and the intersection points. An interface node maybe have two moving
%  direction, just choose a direction which move shorter distance, then the
%  qualities of the elements connecting to this node will decrease less.

interfaceElemIdx= markinterface(elem,phiValue);
interfaceElem = elem(interfaceElemIdx,:);
T = auxstructure(interfaceElem);
edge = double(T.edge);
edge2elem = T.edge2elem;
isSpecial = elemType(interfaceElemIdx(edge2elem(:,1))) & ...
    ~elemType(interfaceElemIdx(edge2elem(:,2))) & edge2elem(:,3) == 1;
isSpecial = isSpecial | (~elemType(interfaceElemIdx(edge2elem(:,1))) & ...
    elemType(interfaceElemIdx(edge2elem(:,2))) & edge2elem(:,4) == 1);
isShortEdge = (edge2elem(:,1) ~= edge2elem(:,2)) & ~(edge2elem(:,3) == 1 & ...
               edge2elem(:,4) == 1) & ~isSpecial;
A = node(edge(isShortEdge,1),:);
B = node(edge(isShortEdge,2),:);
h = sqrt(sum((A-B).^2,2));
M = findintersectbisect(phi,A,B);
h1 = sqrt(sum((A-M).^2,2));
h2 = h - h1;
t = sparse(edge(isShortEdge,[1,2]),edge(isShortEdge,[2,1]),[h1./h,h2./h], size(node,1),size(node,1));
t1 = sparse(edge(isShortEdge,[1,2]),edge(isShortEdge,[2,1]),repmat(cumsum(ones(size(M,1),1)),1,2), size(node,1),size(node,1));
[maxt,I] = max(t,[],1);
isIrregularNode = maxt >= 0.5;
idx = find(isIrregularNode);
idx1 = I(idx)+(idx - 1)*size(node,1);

oldNodeIdx = [oldNodeIdx;idx'];
oldNode = [oldNode;node(idx,:)];

node(idx,:) = M(t1(idx1),:);

phiValue(idx) = 0;

display('Finish Step 4: move the nearest nodes to interface!');
% figure
% showmesh(node,elem);
% findelem(node,elem,interfaceElemIdx,'noindex');
% plotcurve(a);

%% Step 5: Edge swap
%
%  We just consider the case that two elements share a common longest edge,
%  whose one or two of its four vertices are moved on the interface. After
%  moving, the qualities of this two elements will change(and most case,
%  the quality will decrease). If edge swapping can improve the qualities,
%  we do edge swapping for this two elements.
%  
%  But if the common longest edge of the two element is an interface edge,
%  namely its two endpoints are on the interface, we can't use edge
%  swapping. In this case, we deal with it in Step 7. Just find the element
%  patchs of  the first vertices of these two element, and replace this two
%  vertices by the centroids of the element patchs, respectively.

T = auxstructure(elem);
edge = double(T.edge);
edge2elem = T.edge2elem; 
clear T;

leftElemIdx = edge2elem(:,1); % the index of the first neighbor of edge 
rightElemIdx = edge2elem(:,2); % the index of the second neighbor of  edge

isIrregularNode = full(isIrregularNode)';
isLongestEdge = (edge2elem(:,1) ~= edge2elem(:,2)) & ...
    edge2elem(:,3) == 1 & edge2elem(:,4)== 1 & ...
    (isIrregularNode(edge(:,1)) | isIrregularNode(edge(:,2)) | ...
    isIrregularNode(elem(leftElemIdx,1)) | isIrregularNode(elem(rightElemIdx,1)));

isCrossLongestEdge = (edge2elem(:,1) ~= edge2elem(:,2)) & ...
    edge2elem(:,3) == 1 & edge2elem(:,4)== 1 & ...
    isIrregularNode(elem(leftElemIdx,1)) & isIrregularNode(elem(rightElemIdx,1)) &...
    ~isIrregularNode(edge(:,1)) & ~isIrregularNode(edge(:,2));

isInterfaceEdge = isIrregularNode(edge(:,1)) & isIrregularNode(edge(:,2));
isInterfaceEdge = isLongestEdge & isInterfaceEdge ;
isInterfaceEdge = isInterfaceEdge & ~isIrregularNode(elem(leftElemIdx,1)) ;
isInterfaceEdge = isInterfaceEdge & ~isIrregularNode(elem(rightElemIdx,1));

isLongestEdge(isInterfaceEdge) = false;

leftElemIdx = edge2elem(isLongestEdge,1);
rightElemIdx = edge2elem(isLongestEdge,2);

leftElem = elem(leftElemIdx,:);
rightElem = elem(rightElemIdx,:);

leftElemNew = [leftElem(:,2), rightElem(:,1),leftElem(:,1)];
rightElemNew = [rightElem(:,2), leftElem(:,1),rightElem(:,1)];

l1 = triangleratio(node,leftElem);
r1 = triangleratio(node,rightElem);

l2 = triangleratio(node,leftElemNew);
r2 = triangleratio(node,rightElemNew);

notNeedSwap = min([l1,r1],[],2) >= min([l2,r2],[],2) & ~isCrossLongestEdge(isLongestEdge);

leftElemNew(notNeedSwap,:) = leftElem(notNeedSwap,:);
rightElemNew(notNeedSwap,:) = rightElem(notNeedSwap,:);

oldElemIdx =[oldElemIdx;leftElemIdx;rightElemIdx];
oldElem = [oldElem;elem([leftElemIdx;rightElemIdx],:)];

elem(leftElemIdx,:) = leftElemNew;
elem(rightElemIdx,:) = rightElemNew;

display('Finish Step 5: edge swap!')

%% Step 6: Mesh smoothing for near interface nodes
% q = minangle(node,elem);
% [minq,I] = min(q);

T = auxstructure(elem);
% edge = double(T.edge);
neighbor = double(T.neighbor);
clear T;

isInterfaceNode = (msign(phiValue) == 0);

badElem = find(~isInterfaceNode(elem(:,1)) & isInterfaceNode(elem(:,2)) ...
              & isInterfaceNode(elem(:,3)));
v1 = node(elem(badElem,3),:) - node(elem(badElem,2),:);
v2 = node(elem(badElem,1),:) - node(elem(badElem,3),:);
v3 = node(elem(badElem,2),:) - node(elem(badElem,1),:);
badElem = badElem(sum(v2.^2,2)+sum(v3.^2,2) < sum(v1.^2,2)-eps);
moveNodeIdx = elem(badElem,1);
if ~isempty(moveNodeIdx)
   N = size(node,1);
   NT = size(elem,1);
   v2t = sparse(elem,[(1:NT)',(1:NT)',(1:NT)'],1,N,NT);
   area = simplexvolume(node,elem);
   bc = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
   oldNodeIdx = [oldNodeIdx;moveNodeIdx];
   oldNode = [oldNode;node(moveNodeIdx,:)];
   ba = v2t(moveNodeIdx,:)*area;
   node(moveNodeIdx,:) = v2t(moveNodeIdx,:)*([area,area].*bc)./[ba,ba];
end
% case:
badElem2 = find(~isInterfaceNode(elem(:,2)) & isInterfaceNode(elem(:,1)) ...
              & isInterfaceNode(elem(:,3)) & ~elemType(neighbor(:,1)));
v1 = node(elem(badElem2,3),:) - node(elem(badElem2,2),:);
v2 = node(elem(badElem2,1),:) - node(elem(badElem2,3),:);
v3 = node(elem(badElem2,2),:) - node(elem(badElem2,1),:);
badElem2 = badElem2(sum(v2.^2,2)+sum(v3.^2,2)<sum(v1.^2,2)-eps);
% case:
badElem3 = find(~isInterfaceNode(elem(:,3)) & isInterfaceNode(elem(:,1)) ...
              & isInterfaceNode(elem(:,2)) & ~elemType(neighbor(:,1)));
v1 = node(elem(badElem3,3),:) - node(elem(badElem3,2),:);
v2 = node(elem(badElem3,1),:) - node(elem(badElem3,3),:);
v3 = node(elem(badElem3,2),:) - node(elem(badElem3,1),:);
badElem3 = badElem3(sum(v2.^2,2)+sum(v3.^2,2)<sum(v1.^2,2)-eps);
moveNodeIdx = [elem(badElem2,2);elem(badElem3,3)];
if ~isempty(moveNodeIdx)
   if ~exist('v2t','var')
       N = size(node,1);
       NT = size(elem,1);
       v2t = sparse(elem,[(1:NT)',(1:NT)',(1:NT)'],1,N,NT);
   end
   area = simplexvolume(node,elem);
   bc = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
   isMany = sum(v2t(moveNodeIdx,:),2)<=6;
   moveNodeIdx = moveNodeIdx(isMany);
   oldNodeIdx = [oldNodeIdx;moveNodeIdx];
   oldNode = [oldNode;node(moveNodeIdx,:)];
   ba = v2t(moveNodeIdx,:)*area;
   node(moveNodeIdx,:) = v2t(moveNodeIdx,:)*([area,area].*bc)./[ba,ba];  
end
display('Finish Step 6: We got it finally!');
% 
% q = minangle(node,elem);
% [minq,I] = min(q);