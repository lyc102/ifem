node = [1,0; 0,1; -1,0; 0,-1; 0,0];      % nodes
elem = [5,1,2; 5,2,3; 5,3,4];            % elements
bdFlag = [1 0 0; 1 0 0; 1 0 0];
for k = 1:3
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    bdNode = findboundary(elem,bdFlag);
    nodenorm = sqrt(node(bdNode,1).^2+node(bdNode,2).^2);
    node(bdNode,:) = node(bdNode,:)./[nodenorm nodenorm];
    subplot(1,2,1); showmesh(node,elem);                            % plot mesh
    [node,elem] = optmesh(node,elem,1);
    subplot(1,2,2); showmesh(node,elem);                            % plot mesh
end
% findelem(node,elem);                            % plot element indices
% findnode(node,1:5);                             % plot node indices
