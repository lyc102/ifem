
%% Unit Cube 
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
subplot(1,2,1); showmesh3(node,elem);

%% Lshape domain
% mesh
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
[node,elem] = delmesh(node,elem,'x<0 & y<0');
bdFlag = setboundary3(node,elem,'Dirichlet');
% for adaptive mesh refinement, special ordering of elem and HB is needed
[elem,bdFlag,HB] = label3(node,elem,'all',bdFlag);
showmesh3(node,elem);