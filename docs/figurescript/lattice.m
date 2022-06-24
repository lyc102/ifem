function lattice(k) 
[node,elem] = squaremesh([0,1,0,1],1/k);
[node,elem] = delmesh(node,elem,'y>x');
showmesh(node,elem);
findnode(node,'all','noindex','MarkerSize',64);
filename = sprintf('lattice2D_%d.png',k);
saveas(gcf,filename)