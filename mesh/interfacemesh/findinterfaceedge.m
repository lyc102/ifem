function [interiorElem,exteriorElem,interfaceEdge] = findinterfaceedge(node,elem,phi)

NT = size(elem,1);
center = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
isInteriorElem = phi(center)<0;
interiorElem = elem(isInteriorElem,:);
exteriorElem = elem(~isInteriorElem,:);

T = auxstructure(elem);
edge2elem = double(T.edge2elem);

prev = [3;1;2];
next = [2;3;1];

isInterfaceEdge1 = (isInteriorElem(edge2elem(:,1)) & ~isInteriorElem(edge2elem(:,2)));
isInterfaceEdge2 = (isInteriorElem(edge2elem(:,2)) & ~isInteriorElem(edge2elem(:,1)));

interfaceElemIdx = [edge2elem(isInterfaceEdge1,1);edge2elem(isInterfaceEdge2,2)];
localIdx = [edge2elem(isInterfaceEdge1,3);edge2elem(isInterfaceEdge2,4)];

idx1 = (next(localIdx)-1)*NT + interfaceElemIdx;
idx2 = (prev(localIdx)-1)*NT + interfaceElemIdx;
interfaceEdge = [elem(idx1), elem(idx2)];





