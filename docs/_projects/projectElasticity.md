---
permalink: /projects/elasticity/
title: "Project: Finite Element Methods for Linear Elasticity"
sidebar:
    nav: project
---

The objective of this project is to implement finite element methods for solving the linear elasticity equation in two dimensions.

## Part I: Conforming Linear Element

Given a force $\boldsymbol{f}\in L^2(\Omega;\mathbb{R}^3)$, we seek $\boldsymbol{u}\in H_0^1(\Omega; \mathbb{R}^3)$​ such that


$$
2\mu (\nabla ^s \boldsymbol{u}, \nabla ^s \boldsymbol{v}) + \lambda (\text{div} \boldsymbol{u}, \text{div} \boldsymbol{v}) = (\boldsymbol{f}, \boldsymbol{v}) \quad \forall \boldsymbol{v}\in H_0^1(\Omega; \mathbb{R}^3),
$$


where $\lambda$ and $\mu$ are Lamé constants. In materials that are nearly incompressible, the parameter $\lambda \gg 1$ while $\mu = \mathcal{O}(1)$. Utilizing the identity $2 \text{div} \nabla^s \boldsymbol{u} = \Delta \boldsymbol{u} + \text{grad div} \boldsymbol{u}$​, we can derive an equivalent bilinear form


$$
a(\boldsymbol{u}, \boldsymbol{v}) = \mu (\nabla \boldsymbol{u}, \nabla \boldsymbol{v}) + (\lambda + \mu )(\text{div} \boldsymbol{u}, \text{div} \boldsymbol{v}).
$$


We consider $\Omega$ as a unit square. Given a triangulation, the displacement space $\mathbf{u} = (u_1, \; u_2)$ is the linear finite element space.

## Step 1: Assemble the matrix equation

This part involves assembling the stiffness matrix equation, similar to the procedure outlined in  [Project: Linear Finite Element method](https://lyc102.github.io/ifem/projects/fem/).

1. Compute the local stiffness matrix for each element, which is a $6\times 6$ matrix consisting of two copies of $-\Delta$ on the diagonal and the off-diagonal terms $(\text{div} \boldsymbol{\phi}_i, \text{div} \boldsymbol{\phi}_j)$.

   When computing $\rm div\boldsymbol \phi$, expand the linear element into barycentric coordinates and compute the gradient of the basis functions using the function `gradbasis(node,elem)`.

3. Compute the right-hand side of the equation and modify it to incorporate the boundary conditions.

4. Utilize a direct solver to solve the linear algebraic system.

## Step 2: Convergence

Set $\lambda = \mu = 1$ and select a smooth exact solution $\boldsymbol{u}$. Compute the $L^2$ and $H^1$ errors of $\boldsymbol{u} - \boldsymbol{u}_h$.

## Step 3: Locking

For a fixed $h$, choose $\lambda = 10, 100, 1000$ and compute the approximation error. Investigate how the error behaves as $\lambda$ increases to identify any locking phenomena.



## Part II: Locking Free Schemes

In this part, we introduce an artificial pressure $p = \lambda \text{div} \boldsymbol{u}$ and rewrite the equation into perturbed Stokes equations: Find $\boldsymbol u\in H_0^1(\Omega; \mathbb R^3)$ and $p\in L^2_0(\Omega)$​ such that


$$
\begin{aligned}
\mu (\nabla \boldsymbol u, \nabla \boldsymbol v) + (p, {\rm div} \boldsymbol v) & = (\boldsymbol f, \boldsymbol v), &\text{for all } \boldsymbol v\in \boldsymbol H_0^1(\Omega),\\
({\rm div} \boldsymbol u, q) - \frac{1}{\lambda + \mu}(p,q) & = 0, &\text{for all } q\in L_2(\Omega).
\end{aligned}
$$



## Step 1: Stokes pair

Choose either the ${\rm iso}, P_2 - P_0$ or $P_1^{\rm CR} - P_0$ Stokes pair to solve the perturbed Stokes equations. Refer to [Project: Finite Element Methods for Stokes Equations](https://lyc102.github.io/ifem/projects/Stokes/) for details.

## Step 2: Locking free

Verify that the displacement-pressure formulation provides a locking-free approximation. Test the formulation and ensure that there are no locking phenomena present.
