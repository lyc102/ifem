close all;clear
load arctanMesh.mat
figure(1);
subplot(1,2,1)
showmeshpoly(node,elem);
subplot(1,2,2)
showsolutionpoly(node,elem,u);
view(20,60)

set(gcf,'color','w','Position', [100 600 1000 400])
drawnow;