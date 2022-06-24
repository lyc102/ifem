function lattice3(k) 
[node,elem] = cubemesh([0,1,0,1,0,1],1/k);
[node,elem] = delmesh(node,elem,'x>y | z>y | z>x');
showboundary3(node,elem); 
view([10,30]);
axis tight;
findnode3(node,'all','noindex','MarkerSize',64);
filename = sprintf('lattice3D_%d.png',k);
saveas(gcf,filename)