---
permalink: /projects/PBE/
title: "Project: Nonlinear Poisson Boltzmann Equation"
sidebar:
    nav: project
---

The purpose of this project is to implement Newton's method, two-grid method, and nonlinear multigrid method (FAS) for solving the nonlinear elliptic equations. The example is the nonlinear Poisson-Boltzmann equation for the potential $u$​ corresponding to a given charge density

$$
-\Delta u + k^2 \sinh (u) = \rho (x), \quad x\in \Omega\subset \mathbb R^d,
$$
and with Dirichlet boundary condition $u|_{\partial \Omega} = g.$​

We consider a 2-d problem on the unit square $\Omega = (0,1)^2$. Let $\boldsymbol a=(1.0,2.0)/\sqrt{5}$.  We choose $k =1, \rho$, and $g$  such that the exact solution is 

$$
u(\boldsymbol x) = \bar u(0.1+\boldsymbol a\cdot\boldsymbol x),\quad \bar u(s) = \ln \left ( \frac{1+\cos (s)}{1-\cos (s)}\right).
$$


## Step 1: Linearied Poisson Boltzmann Equation

* Given a current approximation $u_k$, derive the linearized Poisson-Boltzmann equation (LPBE) at $u_k$.
* Assemble the matrix equation for the LPBE. Besides the matrix of Laplacian operator, you need to compute the mass matrix corresponding to the L2 inner product. Mass lumping can be applied.

* Use the direct solver \ to solve the matrix equation of LPBE.
* Use a multigrid solver to solve the matrix equation. You can use your own multigrid solver or call `amg` in ifem.



## Step 2: Newton's Method

* Implement the Newton's method. Control the relative error of the residual in the stopping criteria.
* Change the tolerance or maximum iteration steps in the multigrid solver for each LPBE and collect a table of total iteration steps and cpu time for different choices of inner iteration steps.
* Uniform refine the grid several times. List the iteration steps to reach `1e-6` for different $h$​​​ and the approximation error in $H_1(\Omega)$​​​​ norm.



## Step 3: Two-Grid Method

* Apply Newton's method in Step 2 to $H = 1/4$ to obtain a solution $u_H$​.
* Prolongate $u_H$ to a fine space with $h = H^2$  using `nodeinterpolate`.
* Solve one fixed iteration or one Newton's iteration in the fine grid to obtain $u^h$​​
* Change $H$​​ from 1/4 to 1/16 and show the error for $u^h$​​ in $H_1(\Omega)$​​​ norm.
* Compare the computational cost of two-grid method vs Newton's method.



## Step 4: Nonlinear Multigrid FAS

* Implement the nonlinear Gauss-Seidel smoother.

* Test the two level version of FAS.

* Change two level FAS to V-cycle FAS by recrusion.

* Compare the convergence of FAS with Newton's method and Two-Grid method.
