     [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],2);
      u = ones(size(node,1),1);
      [node,elem,~,HB] = bisect3(node,elem,[1 2],[],HB);
      [node,elem,~,HB] = bisect3(node,elem,'all',[],HB);
      [node,elem,~,HB] = bisect3(node,elem,[1 2 3],[],HB);
      [node,elem,~,HB] = bisect3(node,elem,[1 2 3],[],HB);
      showmesh3(node,elem);
      u = nodeinterpolate(u,HB);
      [node,elem,~,HB,indexMap] = coarsen3(node,elem,'all',[],HB);
      u = nodeinterpolate(u,indexMap);
 