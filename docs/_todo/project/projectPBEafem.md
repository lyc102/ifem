# Project: AFEM for Nonlinear Poisson Boltzmann Equation

The purpose of this project is to test Newton's method combined with
adaptivity for solving the nonlinear elliptic equation. The example is
the nonlinear Poisson-Boltzmann equation for the potential u
corresponding to a given charge density $\rho (x)$ reads

$$ -\Delta u + k^2 \sinh (u) = \rho (x) $$

for $x\in \Omega$, and $u|_{\partial \Omega} = g.$

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

* Use the direct solver to solve the matrix equation.

* Use a multigrid solver (e.g. amg) to solve the matrix equation. You can use your own multigrid methods or call `amg` in ifem.

## Step 2: Newton's method on uniform grids

* Implement the Newton's method. Control the relative error of the residual in the stopping criteria.

* Change the tolerance or max iteration steps in multigrid solver and collect a table of total iteration steps and cpu time for different choices of inner iteration.

* Uniform refine the grid and list the iteration steps to reach `1e-6` for different h and compute the approximation error in $H^1$ norm.

## Step 3: Adaptivity

* Run 2D examples: `Lshape`, `crack` and `Kellogg` in iFEM and read the code to learn the usage of AFEM.

* In each Newton iteration, apply the local mesh refinement for the linearized Poisson-Boltzmann equation. You can chose either recovery type or residual type error estimator and the number of local refinement.

* Compare the convergent rate for uniformly refined grids and locally refined grids.

## Step 4: Two-grid Method

* Apply adaptive Newton's method in Step 3 on a coarse grid.
* Refine the grid several times by `uniformrefine`. Get the `HB` output.
* Use `nodeinterpolate` to inerpolate a funciton in coarse grid to the fine grid during the refinement.
* Apply one Newton iteration in the fine grid and possible local refinement.
* Compare the computation time for methods in Step 3 and Step 4.
