---
permalink: /projects/heat/
title: "Project: AFEM for the Heat Equation"
sidebar:
    nav: project
---

The purpose of this project is to implement explict and implicit
numerical methods for solving the parabolic equation. The example is the heat equation 

$$ u_t -\Delta u  = f \quad \text{ with }u |_{\partial \Omega} = g, u(\cdot,0) = u_0.$$

We consider a 2-d problem on the unit square $\Omega = (0,1)^2$ with the
exact solution

$$u(x,t) = \beta (t)\, \exp(-[(x-t+0.5)^2+(y-t+0.5)^2]/0.04),$$ 

with $$\beta (t) = 0.1\,(1-\exp(-10^2(t-0.5)^2)).$$

Adaptive FEM is further applied to capture the singularity of the solution.

## Step 1: Forward Euler, Backward Euler, and Crack-Nicolson methods

- Given a mesh, construct the stiffness matrix `A` for the Laplace operator and the mass matrix `M` for the $L^2$-inner product.

- Given a time step size `dt`, final time `T`, code a for loop over time to involve the solution by either forward, backward Euler or Crack-Nicolson methods.

> Please do not store the approximation at all time steps. Instead only the solution in the previous step `uold` and the current step `u` is needed.

- For implicit methods, use direct solver `A\b` or multigrid solvers to solve the linear system `Au=b`. For meshes generated in `ifem`, `mg(A,b,elem)` is faster than `amg(A,b)`.

## Step 2: Check the Convergence

- Check the convergence rate in time and space. Use the exact solution to get the nodal interpolant `uI` and compute the H1 norm of the error using matrix `A` and the L2 norm using matrix `M`.

- To check the convergence rate in time, fix a small mesh size `h` in space
 and let `dt` vary and vice verse for convergence in space.

## Step 3: Visulization

- Use `showsolution(node,elem,u)` to plot the solution and together with `pause(0.01)` to get an animation. Use `axis` to fix the axis scaling.

> For small time step, do not plot the solution at every time step. Instead plot every 10 or 100 steps.

- You can save the plot into a movie. Read `doc getframe` for an example.

## Step 4: Adaptive Finite Element Methods

- Run 2D examples: `Lshape` in iFEM and read the code to learn the usage of AFEM; see also [Adaptive Finite Element Methods]({{site.baseurl}}{% link _afem/afem.md %})

- In one time step involution, repeat the refinement and coarsen several steps to get a better approximation of the solution. You can control the max iteration steps for AFEM or the maximal number of elements. 

- Use `eta = estimaterecovery(node,elem,u)` instead of residual type error estimators.

- Use `nodeinterpolate` and `eleminterpolate` to interpolate function between different meshes.

- Check the convergence rate for AFEM.

- Make animation for meshes and solutions.


```matlab

```
