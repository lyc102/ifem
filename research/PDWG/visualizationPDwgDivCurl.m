%% plotting result of div curl
% Need a working LaTeX compiler
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

%%
close all;
k = maxIt - 3;

r1 = showrate(h,errL2UTotal,k,'r-d');
r2 = showrate(h,errL2QuTotal,k,'g-d');
r3 = showrate(h,errqlambdaTotal,k,'b-o');
r4 = showrate(h,errSTotal,k,'c-o');

fprintf('rate in e_u:     %6g\n', r1)
fprintf('rate in e_{Qu}:  %6g\n', r2)
fprintf('rate in e_{q,l}: %6g\n', r3)
fprintf('rate in e_s:     %6g\n', r4)

axis tight;
set(gca,'xscale','log','yscale','log');
T = title('Convergence');
XL = xlabel('Mesh size $\log(1/h)$');
YL = ylabel('Error Magnitude');

L = legend('$\Vert {\varepsilon}^{1/2}(\mathbf{u}-\mathbf{u}_h)\Vert$', ...
          ['$h^{' num2str(r1) '}$'],...
          '$\Vert {\varepsilon}^{1/2}(\mathbf{u}_h-Q\mathbf{u})\Vert$', ...
          ['$h^{' num2str(r2) '}$'],...
          '$|\!|\!|(e_{\lambda}, e_{\mathbf{q}}) |\!|\!| $',...
          ['$h^{' num2str(r3) '}$'],...
          '$|\!|\!|e_{s} |\!|\!| $', ...
          ['$h^{' num2str(r4) '}$'],...
       'LOCATION','Best');
set(gcf,'color','w','Position', [200 200 500 500])
set([XL,YL,L],'Interpreter','latex','FontSize', 12);
set(T,'Interpreter','latex','FontSize',16);
grid on;
drawnow;

return;

%% 
close all; %#ok<UNRCH>
figure(1);
option.scale = 3;
option.plot = 'quiver3';
option.stepSize = 0.05;
% option.expr = '~((x>0) & (y<0))'; % L shape
% option.expr = '~(abs(x+0.5)<0.5 & abs(y+0.5)<0.5 & abs(z+0.5)<0.5)'; % cavity
% option.expr = '~(x>-0.5 & x<0 & y<0 & y>-0.5)'; % torus
option.expr = '~((x>-0.5 & x<0 & y<0 & y>-0.5) | (x>0.5 & x<1 & y<0 & y>-0.5))'; % torus 2 holes
option.layer = 4;
center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
    node(elem(:,3),:) + node(elem(:,4),:))/4;
h1 = showvector3(center,pde.exactu(center),option);
% xlabel('x');ylabel('y')

% for cavity
% hold on;
% [nodeInt,elemInt] = cubemesh([-1/2,1/2,-1/2,1/2,-1/2,1/2],1/2);
% nodeInt = nodeInt-1/2;
% showboundary3(nodeInt,elemInt,'x<0', 'facealpha', 0.3, 'edgecolor', 'n', 'edgealpha',0)

axis tight;
axes.SortMethod='ChildOrder';
% zlim([0.1 0.4]); % for L shape example
% view(140,15); % example 1
% view(20,12); % example 2
view(10,15); % example 4
set(gca,'TickLabelInterpreter', 'tex');
drawnow;
set(gcf,'color','w','Position', [100 600 500 500])

%%
figure(2);
option.scale = 2;
option.plot = 'quiver3';
% option.expr = '~((x>0) & (y<0))'; % L shape
% option.expr = '~(abs(x+0.5)<0.5 & abs(y+0.5)<0.5 & abs(z+0.5)<0.5)'; % cavity
% option.expr = '~(x>-0.5 & x<0 & y<0 & y>-0.5)'; % torus
option.expr = '~((x>-0.5 & x<0 & y<0 & y>-0.5) | (x>0.5 & x<1 & y<0 & y>-0.5))'; % torus 2 holes
option.layer = 4;
option.stepSize = 0.05;
center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
    node(elem(:,3),:) + node(elem(:,4),:))/4;
% h2 = showvector3(center,soln.uh2elem,option);
h2 = showvector3(center,harmonic_u,option);
% xlabel('x');ylabel('y')
axis tight;
% zlim([0.1 0.4]); % for L shape example
axes.SortMethod='ChildOrder';
% view(140,15); % example 1
% view(20,12); % example 2
view(10,15); % example 4
set(gcf,'color','w','Position', [100 600 500 500])
set(gca,'TickLabelInterpreter', 'tex');
drawnow;
exportgraphics(gcf,'ex7_uwg.pdf','ContentType','vector')

%% 
figure(3);
errorE2V = accumarray(elem(:),repmat(eqn.errorUelem,4,1),[size(node,1) 1]);
% h = showsolution3(node,elem,errorE2V,'x<0.5','EdgeColor','k'); % example 1
% h = showsolution3(node,elem,errorE2V,'z<0.25 ','EdgeColor','k'); % example 2
h = showsolution3(node,elem,errorE2V,'z>1/4 ','EdgeColor','k'); % example 4
view(-10,55);
zticks([0 0.25])
camproj('perspective')
cb = colorbar; 
caxis([0 0.22])
a=get(cb); %gets properties of colorbar
a = a.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1)+0.05 a(2)+0.03 0.05 0.8]);% To change size
set(gcf,'color','w','Position', [100 600 480 250])