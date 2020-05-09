load meshQuad2Poly.mat
close all;
figure(1);
hold on;
showmeshpoly(node,elem);
findnode(node);
findelem(node,elemQuad);
set(gcf,'Position', [100 1000 800 800])
drawnow;