# Multigrid Solver Unit Test (testMG.m)

## Overview

`testMG.m` is a comprehensive unit test suite for the multigrid solver `mg.m`. It validates the solver's performance across different:

- **Mesh resolutions** (up to 6 levels of uniform refinement)
- **Element types** (P1, P2, P3, CR in 2D; P1, P2 in 3D)
- **Solver variants** (V-cycle, W-cycle, F-cycle, CG, MINRES, GMRES)
- **Convergence rates** (verifies optimal error reduction)

## Usage

### Run All Tests
```matlab
testMG
```

This will execute all 8 test suites and display results in the console.

### Expected Output

The test suite will print detailed tables showing:
- DOF (degrees of freedom) at each refinement level
- Iteration counts for the MG solver
- CPU time comparisons (MG vs Direct solver)
- Error measurements (MG solution vs Direct solution)
- Convergence rates (L2 and H1 norms)

### Test Suite Breakdown

1. **Test 1: P1 Element (2D)** - Linear finite elements
   - Tests 6 refinement levels
   - Generates performance plots
   - Verifies mesh independence

2. **Test 2: P2 Element (2D)** - Quadratic finite elements
   - Tests 5 refinement levels
   - Validates auxiliary space method

3. **Test 3: P3 Element (2D)** - Cubic finite elements
   - Tests 4 refinement levels
   - Validates higher-order elements

4. **Test 4: CR Element (2D)** - Crouzeix-Raviart nonconforming element
   - Tests 5 refinement levels

5. **Test 5: P1 Element (3D)** - Linear tetrahedral elements
   - Tests 4 refinement levels in 3D

6. **Test 6: P2 Element (3D)** - Quadratic tetrahedral elements
   - Tests 3 refinement levels in 3D

7. **Test 7: Different Solvers** - Compares solver variants
   - VCYCLE, WCYCLE, FCYCLE
   - CG, MINRES, GMRES

8. **Test 8: Convergence Rates** - Validates optimal convergence
   - Checks O(h²) for L2 error
   - Checks O(h) for H1 error

## Pass/Fail Criteria

Each test verifies:

✓ **Accuracy**: MG solution matches direct solver (error < 1e-6)  
✓ **Convergence**: Solver reaches tolerance (flag == 0)  
✓ **Mesh Independence**: Iteration count grows sublinearly (< 3x growth)  
✓ **Optimal Rates**: L2 rate ≈ 2.0, H1 rate ≈ 1.0 (for P1 elements)

## Sample Output

```
====================================================
   MULTIGRID SOLVER UNIT TEST
====================================================

----------------------------------------------------
Test 1: P1 Element - 2D Poisson Equation
----------------------------------------------------

Level        DOF        h   MG Iter      MG Time Direct Time      Error
----------------------------------------------------------------------
    1         81   0.2500          7       0.0012     0.0002   8.24e-11
    2        289   0.1250          8       0.0018     0.0008   1.65e-10
    3       1089   0.0625          9       0.0045     0.0135   3.30e-10
    4       4225   0.0312         10       0.0132     0.2845   6.60e-10
    5      16641   0.0156         11       0.0456     5.2341   1.32e-09
    6      66049   0.0078         12       0.1845    89.3256   2.64e-09

Iteration growth factor: 1.33
Test 1 PASSED ✓
```

## Generated Plots

The test suite generates visualization figures:

### Test 1 (P1 Element 2D):
- **CPU Time vs DOF**: Compares MG and Direct solver scaling
- **Iterations vs DOF**: Shows mesh independence
- **Convergence**: Relative residual reduction
- **MG Accuracy**: Error compared to direct solver

### Test 8 (Convergence Rates):
- **Error vs Mesh Size**: Log-log plot showing optimal convergence slopes

## Troubleshooting

### If tests fail:

1. **Path Issues**: Ensure `ifem` is in your MATLAB path
   ```matlab
   addpath(genpath('/path/to/ifem'))
   ```

2. **Missing Dependencies**: Check that required functions exist:
   - `squaremesh`, `cubemesh`
   - `uniformrefine`, `uniformrefine3`
   - `Poisson`, `PoissonP2`, `PoissonP3`, `PoissonCR`
   - `Poisson3`, `Poisson3P2`
   - `sincosdata`, `sincosdata3`

3. **Memory Issues**: For large problems, reduce `maxLevel` in individual test functions

4. **Convergence Failures**: Check `mg.m` options:
   - Increase tolerance: `option.tol = 1e-6`
   - Increase max iterations: `option.maxIt = 200`

## Customization

To modify tests, edit the following in each test function:

```matlab
maxLevel = 5;           % Number of refinement levels
optionMG.tol = 1e-8;    % Convergence tolerance
optionMG.solver = 'CG'; % Solver type
```

## Performance Benchmarks

Typical iteration counts (2D P1 element):

| DOF    | V-cycle | W-cycle | CG w/MG |
|--------|---------|---------|---------|
| 289    | 12      | 8       | 8       |
| 1089   | 14      | 9       | 9       |
| 4225   | 15      | 10      | 10      |
| 16641  | 16      | 11      | 11      |

The iteration count should grow very slowly (logarithmically) with mesh refinement, demonstrating the multigrid method's optimal complexity.

## References

See the main `mg.m` documentation for theoretical background on the multigrid method and auxiliary space preconditioning.
