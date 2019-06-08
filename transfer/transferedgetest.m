%% Check Transfer operator for edge
% The script file to check 2D transfer operator for edges.
% Face element and edge element share the same transfer in 2D.


%% A coarse mesh and several bisection refinement 
node = [0,0; 1,0; 1,1; 0,1];   
elem = [1 2 4; 3 4 2];
[node,elem] = bisect(node,elem,'all');
[node,elem] = bisect(node,elem,[1 2]);
[nodef,elemf,~,~,tree] = bisect(node,elem,[1 6]);

figure(1); subplot(1,2,1)
showmesh(node,elem);
findelem(node,elem);
findnode(node,'all');
figure(1); subplot(1,2,2)
showmesh(nodef,elemf);
findelem(nodef,elemf);
findnode(nodef,'all');

[~,edgef] = dofedge(elemf);
[~,edge] = dofedge(elem);

pro = transferedge(elem,elemf,tree);

u = inline('[4.8 - x(:,2), 5.9 + x(:,1)]','x');
uI_c = edgeinterpolate(u,node,edge);
uI_f = edgeinterpolate(u,nodef,edgef);
u_c2f = pro*uI_c;

% Compare the solution:
disp(' Edges     uI_f        u_c2f      uI_f-u_c2f');
disp(num2str([(1:size(edgef,1))' uI_f u_c2f uI_f-u_c2f]))
disp('the last test is suitable for several possible bisections, so it can used to adaptive mesh.');
disp('if you can see the answer is zero, means our test is successful.');

%% Test for uniform bisect
node = [0,0; 1,0; 1,1; 0,1];   
elem = [1 2 4; 3 4 2];
[nodef,elemf,~,~,tree] = uniformbisect(node,elem);

subplot(1,2,1)
showmesh(node,elem);
findelem(node,elem);
findnode(node,'all');
subplot(1,2,2)
showmesh(nodef,elemf);
findelem(nodef,elemf);
findnode(nodef,'all');

[~,edgef] = dofedge(elemf);
[~,edge] = dofedge(elem);

pro = transferedgebisect(elem,elemf,tree);

u = inline('[4.8 - x(:,2), 5.9 + x(:,1)]','x');
uI_c = edgeinterpolate(u,node,edge);
uI_f = edgeinterpolate(u,nodef,edgef);
u_c2f = pro*uI_c;

% Compare the solution:
disp(' Edges     uI_f        u_c2f      uI_f-u_c2f');
disp(num2str([(1:size(edgef,1))' uI_f u_c2f uI_f-u_c2f]))
disp('the last test is suitable for several possible bisections, so it can used to adaptive mesh.');
disp('if you can see the answer is zero, means our test is successful.');
