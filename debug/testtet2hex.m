%% uniform grid
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
showmesh3(node,elem);
findnode3(node,'all');
[node,hexelem] = tet2hex(node,elem,HB);

%% Adaptive grid
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1/4);
[node,elem] = delmesh(node,elem,'x<0 & y<0');
bdFlag = setboundary3(node,elem,'Neumann');
for k = 1:2
    eta = abs(sign(f(node(elem(:,1),:))) + sign(f(node(elem(:,2),:)))...
            + sign(f(node(elem(:,3),:))) + sign(f(node(elem(:,4),:))));
    refineElem = find(eta < 4);
    [node,elem,bdFlag,HB] = bisect3(node,elem,refineElem,bdFlag,HB);
end
[tempvar,bdFace] = findboundary3(elem,bdFlag);
showmesh(node,bdFace); 
[node,hexelem] = tet2hex(node,elem,HB);

%% sphere
function s = f(p)
    s = sum(p.^2,2) - (0.5)^2;
end