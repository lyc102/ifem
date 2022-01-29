---
permalink: /projects/StokesMAC/
title: "Project: MAC Scheme for Stokes Equations"
sidebar:
    nav: project
---


The purpose of this project is to implement the MAC scheme for solving Stokes equations in two dimensions. 

Reference
* [Progamming of MAC for Stokes Equations](http://math.uci.edu/~chenlong/226/MACcode.pdf)
* [Progamming of Finite Difference Methods](http://math.uci.edu/~chenlong/226/FDMcode.pdf)



## Step 1: Distributive Gauss-Seidel Smoothing (DGS)

1. Gauss-Seidel smoothing of velocity to update `u`;
2. form the residual for the continuity equation `rp = g-Bu`;
3. solve the Poisson equation for pressure `Ap*ep = rp` by G-S;
4. distribute the correction to velocity by `u = u + B'*ep`;
5. update the pressure by `p = p - Ap*ep`.

Every step can be implemented in a matrix-free version.

Then use DGS as an iterative method to solve the Stokes equation for a fixed level. As the level increases, the iteration steps will increase to achieve the same tolerance, e.g. the relative residual norm is less than $10^{-3}$​. 

To debug your code, let the exact solution be zero. Start from a random initial guess. The DGS iteration will push all unknowns to zero. 



## Step 2: Two Level Method

The two level method is

1. Presmoothing by DGS
2. Form residuals for momentum and continunity equations
3. Restrict the residual to the coarse grid
4. Use DGS in the coarse grid till converge (Step 1)
5. Prolongate the correction to the fine grid
6. Postsmoothing by DGS

Note: the index map between coarse and fine grids are slightly different for `u,v,p`.

Test your two level methods for different levels. It should convergence in less than 20 steps and indepedent of the number of levels.



## Step 3: Vcycle Multigrid Method

Recrusively apply the two-level method to the coarse grid problem in the previous step to get a V-cycle method.

* Test the convergence of Vcycle method. Record the iteration steps needed to push the relative residual smaller than a tolerance say $10^{-3}$.

* Compute the error between the computed approximation to the exact solution and show the convergence rate in terms of mesh size `h`. 



## Step 4: Test Examples

1. **Analytic solution.** We use a simple model of colliding flow with analytic solutions to test the code. The domain is $[-1,1]^2$​  Compute the data `f` and Dirichlet boundary condition `g_D` for the analytic solution:
   $$
   u = 20xy^3; \quad v = 5x^4 - 5y^4; \quad p = 60x^2y - 20y^3 + {\rm constant}.
   $$
   Choose the constant s.t. the integral of $p$ over the domain is zero.

2. **Driven cavity problem.** The domain is $[-1,1]^2$​ . Stokes equation with $f=0,g=0$ and zero Dirichlet boundary condition except on the top:

$$
\{ y=1; -1 \leq x \leq 1 \mid u = 1, v = 0 \}.
$$

## 
