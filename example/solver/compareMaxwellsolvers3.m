%% COMPAREMAXWELLSOLVERS3 compare the direct solver with mg and amg
% 
% help getfemmatrix3 fore more options.

%% Cube uniform grids
mesh.shape = 'cube';
mesh.type = 'uniform';
mesh.size = 1e4;
pde = 'Maxwell';
fem = 'ND1';
% get the matrix
[eqn,T] = getfemmatrix3(mesh,pde,fem);
% compare solvers
tic; disp('Direct solver'); 
x1 = eqn.bigA\eqn.f; 
toc;
tic; 
x2 = mgMaxwell(eqn.bigA,eqn.f,eqn.AP,eqn.BP,T.node,T.elem,T.edge,T.HB,...
              eqn.isBdEdge);
toc;
fprintf('Difference between direct and mg solvers %0.2g \n',...
         norm(x1-x2)/norm(eqn.f));

%% Lshape adaptive grids
mesh.shape = 'Lshape';
mesh.type = 'adaptive';
mesh.size = 3e4;
pde = 'Maxwell';
fem = 'ND1';
% get the matrix
[eqn,T,option] = getfemmatrix3(mesh,pde,fem);
% compare solvers
tic; disp('Direct solver'); 
x1 = eqn.bigA\eqn.f; 
toc;
tic; 
option.solver = 'CG';
x2 = mgMaxwell(eqn.bigA,eqn.f,eqn.AP,eqn.BP,T.node,T.elem,T.edge,T.HB,...
              eqn.isBdEdge,option);
toc;
disp(norm(x1-x2));