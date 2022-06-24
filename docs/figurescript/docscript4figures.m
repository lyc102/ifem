[node,elem] = squaremesh([0 1 0 1],0.5);
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[edge, i2, j] = unique(totalEdge,'rows','legacy');
NT = size(elem,1);
elem2edge = reshape(j,NT,3);

%%
showmesh(node,elem);
findedge(node,edge);
findelem(node,elem);
display(elem2edge);
%%
i1(j(3*NT:-1:1)) = 3*NT:-1:1; i1=i1';
k1 = ceil(i1/NT); t1 = i1 - NT*(k1-1);
k2 = ceil(i2/NT); t2 = i2 - NT*(k2-1);
edge2elem = [t1,t2,k1,k2];
%%
showmesh(node,elem);
findelem(node,elem,edge2elem(6,1:2));
findedge(node,edge,6);
display(edge2elem);