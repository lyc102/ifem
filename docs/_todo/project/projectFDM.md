# Project: Finite Difference Methods

The purpose of this project is to implement the finite difference method (5-point stencil) for
solving the Poisson equation in a rectangular domain using matrix-free or tensor product matrix. 

Reference: 

* [Finite Difference Methods](http://math.uci.edu/~chenlong/226/FDM.pdf)
* [Programming of Finite Difference Methods](http://math.uci.edu/~chenlong/226/FDMcode.pdf)

## Step 1: Uniform Grids and Grid Functions

Use `meshgrid` or `ndgrid` to generate a uniform grid of the rectangle $(0,1)\times (0,2)$ with size `h`.

Evaluate a function `f(x,y)` on this uniform grid. Plot it using `mesh` or `surf`.

## Step 2: Matrix-Vector Product

Compare the following three ways of computing `A*u`

* Generate a big sparse matrix using the tensor product of 1-D tri-diagonal finite difference matrix and compute `A*u` using this matrix

* Code matrix-free version of `A*u`

* Use the 1-D matrix and tensor product structure to compute `A*u`

Use them to verify the truncation error. That is, choosing the exact
solution `u` and compute the max norm of `A*u - f`. Change `h` and show the
order of the truncation error.

## Step 3: Boundary Conditions

**Dirichlet boundary condition.** Evaluate the boundary condition at
boundary vertices and move to the right hand side.

**Neumann boundary condition.** Change the stencil near the boundary and
modify the right hand side.

## Step 4: Solve the Linear Algebraic Systems

* Direct methods: Use the big matrix to solve `u = A\f`.

* Iterative methods: Implement Gauss-Seidel method as a subroutine `GS(u,f)` and iterate

      while err > tol
           u = GS(u,f);
           err = norm(f-A*u)/norm(f);
      end

## Step 5: Convergence

* Choose a smooth solution `u` and calculate the right hand side `f` and boundary conditions for the unit square.

* Compute the error `uI - uh` in the maximum norm, where `uI` is the nodal interpolation of `u`.

* Change the mesh size `h` from `1/8,1/16,1/32,1/64` and show the order of the error.
