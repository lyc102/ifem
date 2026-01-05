function testMG
%% TESTMG Unit test for multigrid solver on different mesh resolutions
%
% This function tests the mg.m solver on various:
%   - Mesh resolutions (uniform refinement)
%   - Element types (P1, P2, P3, CR)
%   - Dimensions (2D and 3D)
%   - Different solvers (VCYCLE, WCYCLE, CG, GMRES)
%
% The test verifies:
%   - Convergence to the direct solver solution
%   - Mesh independence (iterations should not grow significantly)
%   - Proper error reduction
%
% Usage:
%   testMG
%
% See also mg, Poisson, Poisson3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clc; close all;

%% Test configurations
fprintf('====================================================\n');
fprintf('   MULTIGRID SOLVER UNIT TEST\n');
fprintf('====================================================\n\n');

% Run all tests
testP1Element2D();
testP2Element2D();
testP3Element2D();
testCRElement2D();
testP1Element3D();
testP2Element3D();
testDifferentSolvers2D();
testConvergenceRates2D();

fprintf('\n====================================================\n');
fprintf('   ALL TESTS COMPLETED SUCCESSFULLY!\n');
fprintf('====================================================\n\n');

end

%% Test 1: P1 Element in 2D
function testP1Element2D()
fprintf('\n----------------------------------------------------\n');
fprintf('Test 1: P1 Element - 2D Poisson Equation\n');
fprintf('----------------------------------------------------\n');

% Initial mesh
[node,elem] = squaremesh([0 1 0 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
pde = sincosdata;

% Test parameters
maxLevel = 6;
N = zeros(maxLevel,1);
h = zeros(maxLevel,1);
iterMG = zeros(maxLevel,1);
errMG = zeros(maxLevel,1);
errDirect = zeros(maxLevel,1);
timeMG = zeros(maxLevel,1);
timeDirect = zeros(maxLevel,1);

fprintf('\n%5s %10s %8s %10s %12s %10s %10s\n', 'Level', 'DOF', 'h', 'MG Iter', 'MG Time', 'Direct Time', 'Error');
fprintf('%s\n', repmat('-',70,1));

for k = 1:maxLevel
    % Refine mesh
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    
    % Assemble system
    option.solver = 'none';
    [soln,eqn] = Poisson(node,elem,bdFlag,pde,option);
    
    % Record mesh info
    N(k) = length(eqn.b);
    h(k) = 1/(sqrt(size(node,1))-1);
    
    % Direct solver (reference solution)
    tic;
    uDirect = eqn.A\eqn.b;
    timeDirect(k) = toc;
    
    % MG solver
    optionMG.solver = 'CG';
    optionMG.tol = 1e-8;
    optionMG.printlevel = 0;
    tic;
    [uMG,info] = mg(eqn.A,eqn.b,elem,optionMG);
    timeMG(k) = toc;
    iterMG(k) = info.itStep;
    
    % Compute error between MG and direct solver
    errMG(k) = norm(uMG-uDirect)/norm(uDirect);
    errDirect(k) = max(info.error(:,2));
    
    % Display results
    fprintf('%5d %10d %8.4f %10d %12.4f %10.4f %10.2e\n', ...
            k, N(k), h(k), iterMG(k), timeMG(k), timeDirect(k), errMG(k));
    
    % Verify convergence
    assert(errMG(k) < 1e-6, 'MG solution differs too much from direct solver!');
    assert(info.flag == 0, 'MG did not converge!');
end

% Check mesh independence (iterations should grow slowly)
if maxLevel > 3
    iterGrowth = iterMG(end) / iterMG(3);
    fprintf('\nIteration growth factor: %.2f\n', iterGrowth);
    assert(iterGrowth < 3, 'Iterations grow too fast with mesh refinement!');
end

% Plot results
figure('Name','P1 Element 2D Test');
subplot(2,2,1);
loglog(N, timeMG, 'r-o', N, timeDirect, 'b-s', 'LineWidth', 2);
xlabel('Number of DOFs'); ylabel('CPU Time (s)');
legend('MG', 'Direct', 'Location', 'northwest');
title('P1 Element: CPU Time vs DOF');
grid on;

subplot(2,2,2);
semilogx(N, iterMG, 'k-o', 'LineWidth', 2);
xlabel('Number of DOFs'); ylabel('MG Iterations');
title('P1 Element: Iterations vs DOF');
grid on;

subplot(2,2,3);
loglog(h, errDirect, 'g-^', 'LineWidth', 2);
xlabel('Mesh size h'); ylabel('Relative Residual');
title('P1 Element: Convergence');
grid on;

subplot(2,2,4);
semilogy(N, errMG, 'm-d', 'LineWidth', 2);
xlabel('Number of DOFs'); ylabel('Error vs Direct Solver');
title('P1 Element: MG Accuracy');
grid on;

fprintf('Test 1 PASSED ✓\n');
end

%% Test 2: P2 Element in 2D
function testP2Element2D()
fprintf('\n----------------------------------------------------\n');
fprintf('Test 2: P2 Element - 2D Poisson Equation\n');
fprintf('----------------------------------------------------\n');

% Initial mesh
[node,elem] = squaremesh([0 1 0 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
pde = sincosdata;

maxLevel = 5;
N = zeros(maxLevel,1);
iterMG = zeros(maxLevel,1);
errMG = zeros(maxLevel,1);

fprintf('\n%5s %10s %10s %12s\n', 'Level', 'DOF', 'MG Iter', 'Error');
fprintf('%s\n', repmat('-',40,1));

for k = 1:maxLevel
    % Refine mesh
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    
    % Assemble system
    option.solver = 'none';
    [soln,eqn] = PoissonP2(node,elem,bdFlag,pde,option);
    
    % Record info
    N(k) = length(eqn.b);
    
    % Direct solver
    uDirect = eqn.A\eqn.b;
    
    % MG solver
    optionMG.solver = 'CG';
    optionMG.tol = 1e-8;
    optionMG.printlevel = 0;
    [uMG,info] = mg(eqn.A,eqn.b,elem,optionMG);
    iterMG(k) = info.itStep;
    
    % Compute error
    errMG(k) = norm(uMG-uDirect)/norm(uDirect);
    
    fprintf('%5d %10d %10d %12.2e\n', k, N(k), iterMG(k), errMG(k));
    
    % Verify
    assert(errMG(k) < 1e-6, 'P2: MG solution differs from direct solver!');
    assert(info.flag == 0, 'P2: MG did not converge!');
end

fprintf('Test 2 PASSED ✓\n');
end

%% Test 3: P3 Element in 2D
function testP3Element2D()
fprintf('\n----------------------------------------------------\n');
fprintf('Test 3: P3 Element - 2D Poisson Equation\n');
fprintf('----------------------------------------------------\n');

% Initial mesh
[node,elem] = squaremesh([0 1 0 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
pde = sincosdata;

maxLevel = 4;
N = zeros(maxLevel,1);
iterMG = zeros(maxLevel,1);
errMG = zeros(maxLevel,1);

fprintf('\n%5s %10s %10s %12s\n', 'Level', 'DOF', 'MG Iter', 'Error');
fprintf('%s\n', repmat('-',40,1));

for k = 1:maxLevel
    % Refine mesh
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    
    % Assemble system
    option.solver = 'none';
    [soln,eqn] = PoissonP3(node,elem,bdFlag,pde,option);
    
    N(k) = length(eqn.b);
    
    % Direct solver
    uDirect = eqn.A\eqn.b;
    
    % MG solver
    optionMG.solver = 'CG';
    optionMG.tol = 1e-8;
    optionMG.printlevel = 0;
    [uMG,info] = mg(eqn.A,eqn.b,elem,optionMG);
    iterMG(k) = info.itStep;
    errMG(k) = norm(uMG-uDirect)/norm(uDirect);
    
    fprintf('%5d %10d %10d %12.2e\n', k, N(k), iterMG(k), errMG(k));
    
    assert(errMG(k) < 1e-6, 'P3: MG solution differs from direct solver!');
    assert(info.flag == 0, 'P3: MG did not converge!');
end

fprintf('Test 3 PASSED ✓\n');
end

%% Test 4: CR Element in 2D
function testCRElement2D()
fprintf('\n----------------------------------------------------\n');
fprintf('Test 4: CR (Crouzeix-Raviart) Element - 2D\n');
fprintf('----------------------------------------------------\n');

% Initial mesh
[node,elem] = squaremesh([0 1 0 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
pde = sincosdata;

maxLevel = 5;
N = zeros(maxLevel,1);
iterMG = zeros(maxLevel,1);
errMG = zeros(maxLevel,1);

fprintf('\n%5s %10s %10s %12s\n', 'Level', 'DOF', 'MG Iter', 'Error');
fprintf('%s\n', repmat('-',40,1));

for k = 1:maxLevel
    % Refine mesh
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    
    % Assemble system
    option.solver = 'none';
    [soln,eqn] = PoissonCR(node,elem,bdFlag,pde,option);
    
    N(k) = length(eqn.b);
    
    % Direct solver
    uDirect = eqn.A\eqn.b;
    
    % MG solver
    optionMG.solver = 'CG';
    optionMG.tol = 1e-8;
    optionMG.printlevel = 0;
    [uMG,info] = mg(eqn.A,eqn.b,elem,optionMG);
    iterMG(k) = info.itStep;
    errMG(k) = norm(uMG-uDirect)/norm(uDirect);
    
    fprintf('%5d %10d %10d %12.2e\n', k, N(k), iterMG(k), errMG(k));
    
    assert(errMG(k) < 1e-6, 'CR: MG solution differs from direct solver!');
    assert(info.flag == 0, 'CR: MG did not converge!');
end

fprintf('Test 4 PASSED ✓\n');
end

%% Test 5: P1 Element in 3D
function testP1Element3D()
fprintf('\n----------------------------------------------------\n');
fprintf('Test 5: P1 Element - 3D Poisson Equation\n');
fprintf('----------------------------------------------------\n');

% Initial mesh
[node,elem] = cubemesh([0 1 0 1 0 1], 0.5);
bdFlag = setboundary3(node,elem,'Dirichlet');
pde = sincosdata3;

maxLevel = 4;
N = zeros(maxLevel,1);
iterMG = zeros(maxLevel,1);
errMG = zeros(maxLevel,1);

fprintf('\n%5s %10s %10s %12s\n', 'Level', 'DOF', 'MG Iter', 'Error');
fprintf('%s\n', repmat('-',40,1));

for k = 1:maxLevel
    % Refine mesh
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    
    % Assemble system
    option.solver = 'none';
    [soln,eqn] = Poisson3(node,elem,bdFlag,pde,option);
    
    N(k) = length(eqn.b);
    
    % Direct solver
    uDirect = eqn.A\eqn.b;
    
    % MG solver
    optionMG.solver = 'CG';
    optionMG.tol = 1e-8;
    optionMG.printlevel = 0;
    [uMG,info] = mg(eqn.A,eqn.b,elem,optionMG);
    iterMG(k) = info.itStep;
    errMG(k) = norm(uMG-uDirect)/norm(uDirect);
    
    fprintf('%5d %10d %10d %12.2e\n', k, N(k), iterMG(k), errMG(k));
    
    assert(errMG(k) < 1e-6, '3D P1: MG solution differs from direct solver!');
    assert(info.flag == 0, '3D P1: MG did not converge!');
end

fprintf('Test 5 PASSED ✓\n');
end

%% Test 6: P2 Element in 3D
function testP2Element3D()
fprintf('\n----------------------------------------------------\n');
fprintf('Test 6: P2 Element - 3D Poisson Equation\n');
fprintf('----------------------------------------------------\n');

% Initial mesh
[node,elem] = cubemesh([0 1 0 1 0 1], 0.5);
bdFlag = setboundary3(node,elem,'Dirichlet');
pde = sincosdata3;

maxLevel = 3;
N = zeros(maxLevel,1);
iterMG = zeros(maxLevel,1);
errMG = zeros(maxLevel,1);

fprintf('\n%5s %10s %10s %12s\n', 'Level', 'DOF', 'MG Iter', 'Error');
fprintf('%s\n', repmat('-',40,1));

for k = 1:maxLevel
    % Refine mesh
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    
    % Assemble system
    option.solver = 'none';
    [soln,eqn] = Poisson3P2(node,elem,bdFlag,pde,option);
    
    N(k) = length(eqn.b);
    
    % Direct solver
    uDirect = eqn.A\eqn.b;
    
    % MG solver
    optionMG.solver = 'CG';
    optionMG.tol = 1e-8;
    optionMG.printlevel = 0;
    [uMG,info] = mg(eqn.A,eqn.b,elem,optionMG);
    iterMG(k) = info.itStep;
    errMG(k) = norm(uMG-uDirect)/norm(uDirect);
    
    fprintf('%5d %10d %10d %12.2e\n', k, N(k), iterMG(k), errMG(k));
    
    assert(errMG(k) < 1e-6, '3D P2: MG solution differs from direct solver!');
    assert(info.flag == 0, '3D P2: MG did not converge!');
end

fprintf('Test 6 PASSED ✓\n');
end

%% Test 7: Different Solver Types
function testDifferentSolvers2D()
fprintf('\n----------------------------------------------------\n');
fprintf('Test 7: Different MG Solver Types (2D P1)\n');
fprintf('----------------------------------------------------\n');

% Setup mesh
[node,elem] = squaremesh([0 1 0 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
pde = sincosdata;

% Refine to reasonable size
for k = 1:4
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
end

% Assemble system
option.solver = 'none';
[soln,eqn] = Poisson(node,elem,bdFlag,pde,option);
uDirect = eqn.A\eqn.b;

% Test different solvers
solvers = {'VCYCLE', 'WCYCLE', 'FCYCLE', 'CG', 'MINRES', 'GMRES'};
fprintf('\n%15s %10s %10s %12s\n', 'Solver', 'Iter', 'Time', 'Error');
fprintf('%s\n', repmat('-',50,1));

for i = 1:length(solvers)
    optionMG.solver = solvers{i};
    optionMG.tol = 1e-8;
    optionMG.printlevel = 0;
    
    tic;
    [uMG,info] = mg(eqn.A,eqn.b,elem,optionMG);
    cpuTime = toc;
    
    errSolver = norm(uMG-uDirect)/norm(uDirect);
    
    fprintf('%15s %10d %10.4f %12.2e\n', ...
            solvers{i}, info.itStep, cpuTime, errSolver);
    
    assert(errSolver < 1e-6, sprintf('%s: Solution error too large!', solvers{i}));
    assert(info.flag == 0, sprintf('%s: Did not converge!', solvers{i}));
end

fprintf('Test 7 PASSED ✓\n');
end

%% Test 8: Convergence Rates
function testConvergenceRates2D()
fprintf('\n----------------------------------------------------\n');
fprintf('Test 8: Error Convergence Rates (2D P1)\n');
fprintf('----------------------------------------------------\n');

% Initial mesh
[node,elem] = squaremesh([0 1 0 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
pde = sincosdata;

maxLevel = 5;
h = zeros(maxLevel,1);
errL2 = zeros(maxLevel,1);
errH1 = zeros(maxLevel,1);

fprintf('\n%5s %10s %12s %12s %10s %10s\n', ...
        'Level', 'h', 'L2 Error', 'H1 Error', 'L2 Rate', 'H1 Rate');
fprintf('%s\n', repmat('-',70,1));

for k = 1:maxLevel
    % Refine mesh
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    h(k) = 1/(sqrt(size(node,1))-1);
    
    % Solve with MG
    option.solver = 'mg';
    option.tol = 1e-10;
    option.printlevel = 0;
    [soln,eqn] = Poisson(node,elem,bdFlag,pde,option);
    
    % Compute errors
    errL2(k) = getL2error(node,elem,pde.exactu,soln.u);
    errH1(k) = getH1error(node,elem,pde.Du,soln.Du);
    
    % Compute rates
    if k > 1
        rateL2 = log(errL2(k-1)/errL2(k))/log(h(k-1)/h(k));
        rateH1 = log(errH1(k-1)/errH1(k))/log(h(k-1)/h(k));
        fprintf('%5d %10.4f %12.2e %12.2e %10.2f %10.2f\n', ...
                k, h(k), errL2(k), errH1(k), rateL2, rateH1);
    else
        fprintf('%5d %10.4f %12.2e %12.2e %10s %10s\n', ...
                k, h(k), errL2(k), errH1(k), '-', '-');
    end
end

% Check optimal convergence rates
% P1 element should have O(h^2) for L2 and O(h) for H1
if maxLevel > 2
    avgRateL2 = mean(log(errL2(2:end-1)./errL2(3:end))./log(h(2:end-1)./h(3:end)));
    avgRateH1 = mean(log(errH1(2:end-1)./errH1(3:end))./log(h(2:end-1)./h(3:end)));
    
    fprintf('\nAverage convergence rates:\n');
    fprintf('  L2: %.2f (expected ~2.0)\n', avgRateL2);
    fprintf('  H1: %.2f (expected ~1.0)\n', avgRateH1);
    
    assert(avgRateL2 > 1.7, 'L2 convergence rate too slow!');
    assert(avgRateH1 > 0.8, 'H1 convergence rate too slow!');
end

% Plot convergence
figure('Name','Convergence Rates');
loglog(h, errL2, 'b-o', h, errH1, 'r-s', ...
       h, h.^2, 'b--', h, h, 'r--', 'LineWidth', 2);
xlabel('Mesh size h');
ylabel('Error');
legend('L2 error', 'H1 error', 'O(h^2)', 'O(h)', 'Location', 'southeast');
title('Convergence Rates for P1 Element');
grid on;

fprintf('Test 8 PASSED ✓\n');
end
