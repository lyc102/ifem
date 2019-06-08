%% Check Transfer operator for edge
% The script file to check 2D transfer operator for edge.

%% Test adaptive grids
clear all;
node = [0,0; 1,0; 1,1; 0,1];   
elem = [2,3,1; 4,1,3];
for i = 1:2
    [node,elem] = bisect(node,elem,'all');
end
[nodef,elemf] = bisect(node,elem,[1 4 5]);
[nodec,elemc,~,~,tree] = coarsen(nodef,elemf,'all');
[elem2edgec,edgec] = dofedge(elemc);
[elem2edgef,edgef] = dofedge(elemf);

figure(1)
showmesh(nodec,elemc);
% findnode(nodec);
findelem(nodec,elemc);
findedge(nodec,edgec);

figure(2)
showmesh(nodef,elemf);
findedge(nodef,edgef);
findelem(nodef,elemf);

pro = transferedgecoarsen(elemc,elemf,tree);

u = inline('[3 - x(:,2), 4 + x(:,1)]','x');
uI_c = edgeinterpolate(u,nodec,edgec);
uI_f = edgeinterpolate(u,nodef,edgef);

u_c2f = pro*uI_c;

disp('begin test :')
disp(uI_f - u_c2f);

%% Test boundary edges
bdFlag = setboundary(nodef,elemf,'Dirichlet');
isFreeEdge(elem2edgef(bdFlag == 0)) = true;
freeEdge = find(isFreeEdge);

[pro,freeEdgec] = transferedgecoarsen(elemc,elemf,tree,freeEdge);

u = inline('[sin(pi*x(:,1)).*sin(pi*x(:,2)), sin(pi*x(:,1)).*sin(pi*x(:,2))]','x');
uI_c = edgeinterpolate(u,nodec,edgec);
uI_f = edgeinterpolate(u,nodef,edgef);
u_c2f = zeros(size(edgef,1),1);
u_c2f(freeEdge) = pro*uI_c(freeEdgec);

disp('begin test :')
disp(uI_f - u_c2f);