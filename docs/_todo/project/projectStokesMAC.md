# Project: MAC Scheme for Stokes Equations

The purpose of this project is to implement the simple and popular MAC scheme for solving Stokes equations in two dimensions. 

Reference
* [Progamming of MAC for Stokes equations](http://math.uci.edu/~chenlong/226/MACcode.pdf)

## Test Examples

* Analytic solution. We use a simple model of colliding flow with
analytic solutions to test the code. The domain is $[-1,1]^2$  Compute the data `f` and Dirichlet boundary condition `g_D` for the analytic solution:

$$u = 20xy^3; \quad v = 5x^4 - 5y^4; \quad p = 60x^2y - 20y^3 + {\rm constant}.$$

* Driven cavity problem. The domain is [-1,1]^2. Stokes equation with zero Dirichlet boundary condition except on the top:

$$ \{ y=1; -1 <= x <= 1 | u = 1, v = 0 \}.$$

## Step 1: Gauss-Seidel smoothing of velocity

Given a pressure approximation, relax the momentum equation to update velocity. See [Project: Multigrid Methods](projectMG.html) on the matrix free implemenation of G-S relaxation.

Note that the boundary or near boundary dof should be updated diffeently. The stencil should be changed according to different boundary conditions.

# Step 2: Distributive Gauss-Seidel Smoothing (DGS)

1. Gauss-Seidel smoothing of velocity to update `u`;
* form the residual for the continuity equation `rp = g-Bu`;
* solve the Poisson equation for pressure `Ap*ep = rp` by G-S;
* distribute the correction to velocity by `u = u + B'*ep`;
* update the pressure by `p = p - Ap*ep`.

Every step can be implemented in a matrix-free version; see  [Progamming of MAC for Stokes equations](http://math.uci.edu/~chenlong/226/MACcode.pdf).

Use DGS as an iterative method to solve the Stokes equation. The iteration steps could be very large but will converges.

# Step 3: Two level method

The two level method is

1. presmoothing by DGS
* form residuals for momentum and continunity equations
* restrict the residual to the coarse grid
* iterate DGS in the coarse grid till converge
* prolongate the correction to the fine grid
* postsmoothing by DGS

Note: the index map between coarse and fine grids are slightly different for `u,v,p`; see  [Progamming of MAC for Stokes equations](http://math.uci.edu/~chenlong/226/MACcode.pdf).

Test your two level methods for different levels. It should convergence in less than 20 steps and indepedent of the number of levels.

## Step 4: Vcycle multigrid method

Recrusively apply the two-level method to the coarse grid problem
in the previous step to get a V-cycle method.

* Test the convergence of Vcycle method. Record the iteration steps needed to push the relative residual smaller than a tolerance.

* Compute the error between the computed approximation to the exact solution and show the convergence rate in terms of mesh size `h`. 
