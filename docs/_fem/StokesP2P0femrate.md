# P2P0 Finite Elements for Stokes Equations

This example is to show the convergence of P2-P0 finite elements for the Stokes equation on the unit square:

$$- \Delta u + {\rm grad}\, p  = f \quad {\rm div}\, u    = 0  \quad  \text{ in } \quad \Omega,$$

with the pure Dirichlet boundary condition. The solver is based on a DGS type smoother. 

**References**:
- [Finite Element Methods For Stokes Equations](http://www.math.uci.edu/~chenlong/226/FEMStokes.pdf)
- [Project: Finite Element Methods for Stokes Equations](../project/projectFEM.html)

**Subroutines**:

    - StokesP2P0
    - squareStokes
    - femStokes
    - Stokesfemrate
    
The method is implemented in `StokesP2P0` subroutine and can be tested in `squareStokes`. Together with other elements (P2P0, P2P1, isoP2P0, isoP2P1, P1bP1), `femStokes` provides a concise interface to solve Stokes equation. The P2-P0 element is tested in `Stokesfemrate`. 

## P2-P0 element

The velocity is P2 Lagrange element and the pressure is P0 piecewise constant element. 

We plot the dof below and refer to [PoissonP2femrate](PoissonP2femrate.html) for basis and data structure for P2 element on triangles.


```matlab
clear all;
imatlab_export_fig('print-png')  % Static png figures.
%% Local indexing of DOFs
node = [0,0; 1,0; 0.5, sqrt(3)/2];
elem = [1 2 3];
edge = [2 3; 1 3; 1 2];
% elem2dof = 1:6;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.3,0.2]);
subplot(1,2,1)
showmesh(node,elem);
findnode(node);
findedgedof(node,edge);
subplot(1,2,2)
showmesh(node,elem);
findelem(node,elem);
```


    
![png](StokesP2P0femrate_files/StokesP2P0femrate_3_0.png)
    


## Dirichlet boundary condition


```matlab
%% Setting
% mesh
[node,elem] = squaremesh([0,1,0,1],0.25);
[node,elem] = uniformrefine(node,elem);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% pde
pde = Stokesdata1; 
% options
option.L0 = 0;
option.maxIt = 4;
option.printlevel = 1;
option.solver = 'mg';
```


```matlab
option.elemType = 'P2P0';
femStokes(mesh,pde,option);
```

    #dof:   2434,  #nnz:  25570, level:  3  MG WCYCLE iter: 12,  err = 2.2087e-09,  time = 0.16 s
    #dof:   9986,  #nnz: 108386, level:  4  MG WCYCLE iter: 12,  err = 3.3394e-09,  time = 0.25 s
    #dof:  40450,  #nnz: 446050, level:  5  MG WCYCLE iter: 12,  err = 3.8477e-09,  time =  1.1 s
    #dof: 162818,  #nnz: 1809506, level:  6  MG WCYCLE iter: 12,  err = 4.0427e-09,  time =  5.3 s
    Table: Error
     #Dof        h      |u_I-u_h|_1    ||u-u_h||   ||u_I-u_h||_{max}
    
      2690   6.25e-02   6.64680e-01   1.21603e-02   3.81899e-02
     10498   3.12e-02   3.40775e-01   3.17639e-03   1.02211e-02
     41474   1.56e-02   1.72404e-01   8.10961e-04   2.67100e-03
    164866   7.81e-03   8.66869e-02   2.04806e-04   6.91748e-04
    
     #Dof        h      ||p_I-p_h||    ||p-p_h||   
    
      2690   6.25e-02   1.08871e-01   7.05209e-01
     10498   3.12e-02   3.49510e-02   3.50301e-01
     41474   1.56e-02   1.06474e-02   1.74622e-01
    164866   7.81e-03   3.15916e-03   8.72085e-02
    
    Table: CPU time
     #Dof    Assemble     Solve      Error      Mesh    
    
      2690   8.00e-02   1.60e-01   4.00e-02   1.00e-02
     10498   3.00e-02   2.50e-01   2.00e-02   0.00e+00
     41474   1.00e-01   1.07e+00   8.00e-02   1.00e-02
    164866   4.40e-01   5.33e+00   1.50e-01   2.00e-02
    



    
![png](StokesP2P0femrate_files/StokesP2P0femrate_6_1.png)
    



    
![png](StokesP2P0femrate_files/StokesP2P0femrate_6_2.png)
    


## Conclusion

Optimal order convergence of velocity and pressure is observed. First order for velocity in H1 norm and for pressure in L2 norm. Second order for velocity in L2 and maximum norm. For pressue, superconvergence between interpolate pI and ph is observed.

Multigrid solver based on DGS smoother converges uniformly.

To-Do: test other boundary conditions.
