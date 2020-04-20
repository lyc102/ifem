## coarsen3

- In the input and output, the ordering of arguments is changed to `bdFlag, HB` 
- To use **coarsen3** , `HB` is needed. Therefore add setup in **label3**  to initialize `HB`. 



## Version Changes

Functions in differenti versions of MATLAB might be different. To be compatible, use version check for `unique` function and `DelaunayTri` (earlier version) v.s. `delaunayTriangulation` (2014 and after).

Search  `matlabversion` in iFEM folder to find all functions affected. 



## unique

In matlab, the behavior of `unique`  has been changed.  This includes:

      -	occurrence of indices in IA and IC switched from last to first
      -	IA and IC will always be column index vectors
    
    If this change in behavior has adversely affected your code, you may 
    preserve the previous behavior with:
      
         [C,IA,IC] = unique(A,'legacy')

In ifem, we need 'last'. So in `myunique` function different call of `unique` is provided depending on the version of matlab. 



## Neumann problem

Previously I pitch one point and set `u(1) = 0` , i.e. `freeNode = 2:N` and `fixedNode =1`. After that, enforce $\int_{\Omega} u = 0$ by a constant shift.  

But this requires the exact solution satisfies the condition $u(x(1))=0$ which then fails for some data. In 2D, `Poissonfemrate` use `sincosdata` will give rate 1.8 while change to `sincosNeuman` will give optimal order 2. 

In 3-D, the problem is worse. Rate 1.5 for Neuman and 1.3 for Robin for the error $\|u_I - u_h\|_{\infty}$ . 

So in the new formulation, I perturbe the matrix by `A(1,1) = A(1,1) + 1e-6`. Then the matrix is non-singular. But the same problem still exists. Even change to a consistent data. 

One reason I suspect is the ill-conditioning of the sub-matrix. In 3-D, it is $h^{-3}$ while in 2D is $h^{-2}(1+|\log h|)$.

 **Reference**

- P.B. Bochev, R.B. Lehoucq, On the Finite Element Solution of the Pure Neumann Problem, SIAM Rev. 47 (2005) 50–66.

My MG solver works pretty well but the ill-condition brings an issue to the perturbated problem.

The problem is still the same. We didn't find a solution to $Au = b$ but a perturbated one. The perturbation is small in average (L2 type norm) but can be a problem for maximum problem. 

Here is a simple perturbed analysis. We are solving $ (A + \epsilon) u_{\epsilon} = b$. Substract the equation $(A+\epsilon )u = b + \epsilon u$ to get the equation $ (A + \epsilon) e = \epsilon u$ and thus $\|e\|\leq \| (A+\epsilon)^{-1}\|\|\epsilon u\|$. The matrix $A$ is singular and if $\epsilon < \lambda_2$,  the norm $\| (A+\epsilon)^{-1}\| = \epsilon^{-1}$ and no control of the error. This is the stability issue due to the ill-conditioned matrix. 

To solve the problem, we need to formulate a saddle point system to enforce the constaint $\int_{\Omega}u =0$ into the system which is not straight forward. 

**Conclusion** Accept the current treatment and be aware that for pure Neumann problem, the rate of the computed solution in maximum norm could be degenerated slightly especially in 3D. 

