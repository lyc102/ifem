
%% Packman: refine and project boundary to the circle
node = [1,0; 0,1; -1,0; 0,-1; 0,0];      % nodes
elem = [5,1,2; 5,2,3; 5,3,4];            % elements
bdFlag = [1 0 0; 1 0 0; 1 0 0];
for k = 1:4
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    bdNode = findboundary(elem,bdFlag);
    nodenorm = sqrt(node(bdNode,1).^2+node(bdNode,2).^2);
    node(bdNode,:) = node(bdNode,:)./[nodenorm nodenorm];
    [node,elem] = optmesh(node,elem,1);
end
figure; showmesh(node,elem);

%% Packman: unstructured mesh
fd = inline('ddiff(dcircle(p,0,0,1),drectangle(p,0,1,-1,0))','p');
box=[-1,-1;1,1];
fix=[0,-1; 0,0; 1,0];
[node,elem]=odtmesh2d(fd,@huniform,0.05,box,fix,1);
figure; showmesh(node,elem);

%% Lshape domain
[node,elem] = squaremesh([-1,1,-1,1],0.5);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
showmesh(node,elem);

%% Crack domain
node = [1,0; 0,1; -1,0; 0,-1; 0,0; 1,0];        % nodes
elem = [5,1,2; 5,2,3; 5,3,4; 5,4,6];            % elements
[node,elem] = uniformbisect(node,elem);
showmesh(node,elem);                            % plot mesh
line([node(5,1); node(6,1)'],[node(5,2)'; node(6,2)'],'LineWidth',2,'Color','r');

%% Circle mesh
[node,elem] = circlemesh(0,0,1,0.4);
showmesh(node,elem);
[node,elem] = uniformrefine(node,elem);
[node,elem] = uniformrefine(node,elem);
bdNode = findboundary(elem);
nodenorm = sqrt(node(bdNode,1).^2+node(bdNode,2).^2);
node(bdNode,:) = node(bdNode,:)./[nodenorm nodenorm];
[node,elem] = optmesh(node,elem,1);
showmesh(node,elem);