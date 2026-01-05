%% Quick Test Example for mg.m
% This script demonstrates basic usage of the MG solver
% and can be used for quick validation

clear; clc; close all;

fprintf('====================================================\n');
fprintf('   Quick MG Solver Test\n');
fprintf('====================================================\n\n');

%% Example 1: Simple P1 Element Test
fprintf('Example 1: P1 Element - Compare MG vs Direct Solver\n');
fprintf('----------------------------------------------------\n');

% Create mesh
[node,elem] = squaremesh([0 1 0 1],0.5);
for k = 1:5
    [node,elem] = uniformrefine(node,elem);
end

% Setup PDE
pde.f = 1;
pde.g_D = 0;
bdFlag = setboundary(node,elem,'Dirichlet');

% Assemble system
option.solver = 'none';
[soln,eqn] = Poisson(node,elem,bdFlag,pde,option);

fprintf('Number of unknowns: %8.0u\n',length(eqn.b))
fprintf('Matrix size: %8.0u x %8.0u\n', size(eqn.A,1), size(eqn.A,2))
fprintf('Number of nonzeros: %8.0u\n\n', nnz(eqn.A))

% Direct solver
tic; 
fprintf('Running direct solver...\n');
x1 = eqn.A\eqn.b; 
t_direct = toc;
fprintf('  Time: %.4f s\n\n', t_direct);

% MG solver
tic; 
fprintf('Running multigrid solver...\n');
optionMG.solver = 'CG';
optionMG.printlevel = 1;
x2 = mg(eqn.A,eqn.b,elem,optionMG); 
t_mg = toc;

% Compare results
fprintf('\nSpeedup: %.2fx\n', t_direct/t_mg);
fprintf('Difference between direct and mg solvers: %0.2e\n',...
         norm(x1-x2)/norm(eqn.b));

%% Example 2: Test different solvers
fprintf('\n\n====================================================\n');
fprintf('Example 2: Compare Different MG Solver Types\n');
fprintf('====================================================\n\n');

solverTypes = {'VCYCLE', 'WCYCLE', 'CG', 'GMRES'};
fprintf('%10s | %8s | %10s | %12s\n', 'Solver', 'Iter', 'Time (s)', 'Error');
fprintf('%s\n', repmat('-', 45, 1));

for i = 1:length(solverTypes)
    optionMG.solver = solverTypes{i};
    optionMG.printlevel = 0;
    optionMG.tol = 1e-8;
    
    tic;
    [x_mg, info] = mg(eqn.A, eqn.b, elem, optionMG);
    t = toc;
    
    err = norm(x1 - x_mg)/norm(x1);
    
    fprintf('%10s | %8d | %10.4f | %12.2e\n', ...
            solverTypes{i}, info.itStep, t, err);
end

%% Example 3: Test on different mesh resolutions
fprintf('\n\n====================================================\n');
fprintf('Example 3: Mesh Independence Test\n');
fprintf('====================================================\n\n');

fprintf('%8s | %10s | %8s | %10s\n', 'DOF', 'MG Iter', 'Time', 'Residual');
fprintf('%s\n', repmat('-', 40, 1));

[node,elem] = squaremesh([0 1 0 1],0.5);
for k = 1:6
    [node,elem] = uniformrefine(node,elem);
    
    option.solver = 'none';
    [~,eqn] = Poisson(node,elem,[],pde,option);
    
    optionMG.solver = 'CG';
    optionMG.printlevel = 0;
    optionMG.tol = 1e-8;
    
    tic;
    [x_mg, info] = mg(eqn.A, eqn.b, elem, optionMG);
    t = toc;
    
    fprintf('%8d | %10d | %8.4f | %10.2e\n', ...
            length(eqn.b), info.itStep, t, info.stopErr);
end

%% Example 4: 3D Problem
fprintf('\n\n====================================================\n');
fprintf('Example 4: 3D Poisson Problem\n');
fprintf('====================================================\n\n');

% Create 3D mesh
[node,elem] = cubemesh([0 1 0 1 0 1], 0.5);
for k = 1:3
    [node,elem] = uniformrefine3(node,elem);
end

% Setup 3D PDE
pde3d = sincosdata3;
bdFlag = setboundary3(node,elem,'Dirichlet');

% Assemble
option.solver = 'none';
[soln,eqn] = Poisson3(node,elem,bdFlag,pde3d,option);

fprintf('3D Problem size: %8.0u DOFs\n', length(eqn.b));

% Solve with MG
optionMG.solver = 'CG';
optionMG.printlevel = 1;
tic;
x3d = mg(eqn.A, eqn.b, elem, optionMG);
t3d = toc;

fprintf('\nTotal solution time: %.4f s\n', t3d);

fprintf('\n====================================================\n');
fprintf('   All examples completed!\n');
fprintf('====================================================\n\n');
