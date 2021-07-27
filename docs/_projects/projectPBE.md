---
permalink: /projects/PBE/
title: "Project: Nonlinear Poisson Boltzmann Equation"
sidebar:
    nav: project
---

The purpose of this project is to implement Newton's method, two-grid and nonlinear multigrid method (FAS) for solving the nonlinear elliptic equation. The example is the nonlinear
Poisson-Boltzmann equation for the potential $u$ corresponding to a given
charge density reads

$$ -\Delta u + k^2 \sinh (u) = \rho (x) $$

for $x\in \Omega\subset \mathbb R^d$, and with Dirichlet boundary condition $u|_{\partial \Omega} = g.$

For $k = 1$ and $\rho = 0$, an exact solution in 1-d is given by 

$$ \bar u(s) = \ln \left ( \frac{1+\cos (s)}{1-\cos (s)}\right).$$ 

We consider a 2-d problem on the unit square $\Omega = (0,1)^2$. Let
$\boldsymbol a=(1.0,2.0)/\sqrt{5}$. We choose $k =1, \rho$, and $g$ such that
the exact solution is 

$$u(\boldsymbol x) = \bar u(0.1+\boldsymbol a\cdot\boldsymbol x).$$

## Step 1: Linearied Poisson Boltzmann Equation

* Given a current approximation of u, derive the linearized Poisson-Boltzmann equation (LPBE) at u.

* Assemble the matrix equation for the LPBE. Besides the matrix of Laplacian operator, you need to compute the mass matrix corresponding to the L2 inner product. You can use three vertices quadrature rule i.e.

    $$\int _{\tau} f(x) dx = \frac{1}{3}\sum _{i=1}^3f(x_i)|\tau|.$$ 

    Then the mass matrix becomes diagonal. This is known as mass lumping.

* Use the direct solver \ to solve the matrix equation.

* Use a multigrid solver (e.g. amg) to solve the matrix equation. You can use your own multigrid solver or call `amg` in ifem.

## Step 2: Newton's Method

* Implement the Newton's method. Control the relative error of the residual in the stopping criteria.

* Change the tolerance or maximum iteration steps in the multigrid solver and collect a table of total iteration steps and cpu time for different choices of inner iteration steps.

* Uniform refine the grid and list the iteration steps to reach `1e-6` for different h and compute the approximation error in $H_1(\Omega)$ norm.

## Step 3: Two-Grid Method

* Apply Newton's method in Step 2 to $H = 1/4$ to obtain a solution $u_H$

* Prolongate $u_H$ to a fine space with $h = H^2$ by using the prolongation matrix or subroutine [Project: Multigrid Methods](projectMG.html)

* Solve one fixed iteration or one Newton's iteration in the fine grid to obtain $u_h$

* Change $H$ from 1/4 to 1/16 and show the error for $u_h$ in $H_1(\Omega)$ norm.

## Step 4: Nonlinear Multigrid: FAS

* Implement the nonlinear Gauss-Seidel smoother.

* Test the two level version of FAS.

* Change two level FAS to V-cycle FAS by recrusion.

* Compare the convergence of FAS with Newton's method and Two-Grid method.
