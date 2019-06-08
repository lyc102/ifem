function showmatrix(A,node,theta)
%% SHOWMATRIX displays a matrix.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('node','var')
   [node,D] = eigs(A,2,'SM'); 
end
if ~exist('theta','var'), theta = 0; end
hold on
[i,j,aij] = find(A);
edge = [i j];
idx = (abs(aij) >= theta);
figure;
showgraph(node,edge(idx,:));
% line([node(i,1)'; node(j,1)'],[node(i,2)'; node(j,2)'],'Color','k');
nnz = length(i);
u = zeros(nnz,1);
v = u;
x = 0.5*node(i,1)+0.5*node(j,1); 
y = 0.5*node(i,2)+0.5*node(j,2);
if (size(node,2)==3) % 3-D
    z = 0.5*node(i,3)+0.5*node(j,3);
else   % 2-D
    z = zeros(nnz,1);
end
idx = ((aij<0) & (abs(aij)>theta*mean(abs(aij))));
plot(x(idx), y(idx), 'g*', 'MarkerSize', 6);
idx = ((aij>0) & (abs(aij)>theta*mean(abs(aij))));
plot(x(idx), y(idx), 'c*', 'MarkerSize', 6);
quiver3(x,y,z,u,v,aij,'lineWidth',2);
%text(x+0.002,y+0.002,z+0.01,num2str(aij,3),'FontSize',14);
hold off