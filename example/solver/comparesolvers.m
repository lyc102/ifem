%% COMPARESOLVERS compare the direct solver with mg and amg
% 
% help getfemmatrix fore more options.

%% Square uniform grids
mesh.shape = 'square';
mesh.type = 'uniform';
mesh.size = 2e5;
pde = 'Poisson';
fem = 'P2';
% get the matrix
[eqn,T] = getfemmatrix(mesh,pde,fem);
% compare solvers
tic; disp('Direct solver'); x1 = eqn.A\eqn.b; toc;
tic; x2 = mg(eqn.A,eqn.b,T.elem); toc;
tic; x3 = amg(eqn.A,eqn.b); toc;
fprintf('Difference between direct and mg, amg solvers %0.2g, %0.2g \n',...
         norm(x1-x2)/norm(eqn.b),norm(x1-x3)/norm(eqn.b));
     
%% Lshape adaptive grids
mesh.shape = 'Lshape';
mesh.type = 'adaptive';
mesh.size = 4e4;
pde = 'Poisson';
fem = 'P1';
% get the matrix
[eqn,T] = getfemmatrix(mesh,pde,fem);
% compare solvers
tic; disp('Direct solver'); x1 = eqn.A\eqn.b; toc;
tic; x2 = mg(eqn.A,eqn.b,T.elem); toc;
tic; x3 = amg(eqn.A,eqn.b); toc;
fprintf('Difference between direct and mg, amg solvers %0.2g, %0.2g \n',...
         norm(x1-x2)/norm(eqn.b),norm(x1-x3)/norm(eqn.b));