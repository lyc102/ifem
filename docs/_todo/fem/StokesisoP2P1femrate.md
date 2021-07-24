# isoP2P1 Finite Elements for Stokes Equations

This example is to show the convergence of isoP2-P1 finite elements for the Stokes equation on the unit square:

$$- \Delta u + {\rm grad}\, p  = f \quad {\rm div}\, u    = 0  \quad  \text{ in } \quad \Omega,$$

with the pure Dirichlet boundary condition. The solver is based on a DGS type smoother. 

**References**:
- [Finite Element Methods For Stokes Equations](http://www.math.uci.edu/~chenlong/226/FEMStokes.pdf)
- [Project: Finite Element Methods for Stokes Equations](../project/projectFEM.html)

**Subroutines**:

    - StokesisoP2P1
    - squareStokes
    - femStokes
    - Stokesfemrate
    
The method is implemented in `StokesisoP2P1` subroutine and can be tested in `squareStokes`. Together with other elements (P2P0, P2P1, isoP2P0, isoP2P1, P1bP1), `femStokes` provides a concise interface to solve Stokes equation. The isoP2-P1 element is tested in `Stokesfemrate`.

## isoP2-P1 element

The velocity is P1 Lagrange element but on a uniform refined mesh and the pressure is P1 Lagrange element. Therefore standard linear elements can be used. 


```matlab
imatlab_export_fig('print-png')  % Static png figures.
node = [0,0; 1,0; 0.5, sqrt(3)/2];
elem = [1, 2, 3];
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.3,0.2]);
subplot(1,2,2); showmesh(node,elem); findnode(node,'all','noindex');
[node,elem] = uniformrefine(node,elem);
subplot(1,2,1); showmesh(node,elem); findnode(node,'all','noindex');
```


    
![png](StokesisoP2P1femrate_files/StokesisoP2P1femrate_3_0.png)
    


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
option.elemType = 'isoP2P1';
femStokes(mesh,pde,option);
```

    #dof:   2211,  #nnz:  24940, level:  3  MG WCYCLE iter:  7,  err = 8.1294e-10,  time = 0.11 s
    #dof:   9027,  #nnz: 103788, level:  4  MG WCYCLE iter:  7,  err = 1.4437e-09,  time = 0.15 s
    #dof:  36483,  #nnz: 423276, level:  5  MG WCYCLE iter:  7,  err = 1.4042e-09,  time = 0.47 s
    #dof: 146691,  #nnz: 1709420, level:  6  MG WCYCLE iter:  7,  err = 1.1486e-09,  time =  2.1 s
    Table: Error
     #Dof        h      |u_I-u_h|_1    ||u-u_h||   ||u_I-u_h||_{max}
    
      2467   6.25e-02   6.44313e-02   3.76942e-03   1.82827e-02
      9539   3.12e-02   2.02371e-02   9.25164e-04   4.71869e-03
     37507   1.56e-02   6.56641e-03   2.29540e-04   1.19965e-03
    148739   7.81e-03   2.19661e-03   5.71966e-05   3.02509e-04
    
     #Dof        h      ||p_I-p_h||    ||p-p_h||   
    
      2467   6.25e-02   2.20927e-01   1.55529e-01
      9539   3.12e-02   6.03054e-02   4.26261e-02
     37507   1.56e-02   1.67145e-02   1.18215e-02
    148739   7.81e-03   4.80784e-03   3.38912e-03
    
    Table: CPU time
     #Dof    Assemble     Solve      Error      Mesh    
    
      2467   1.30e-01   1.10e-01   5.00e-02   1.00e-02
      9539   4.00e-02   1.50e-01   2.00e-02   0.00e+00
     37507   1.50e-01   4.70e-01   8.00e-02   0.00e+00
    148739   5.20e-01   2.06e+00   1.70e-01   1.00e-02
    



    
![png](StokesisoP2P1femrate_files/StokesisoP2P1femrate_6_1.png)
    



    
![png](StokesisoP2P1femrate_files/StokesisoP2P1femrate_6_2.png)
    


## Conclusion

Convergence order for velocity in L2 and maximum norm is optimal (2nd order). Half order superconvergence between nodal interpolate uI and uh is observed. As pressure is approximated by linear element, it is almost second order although theoretically it is only first order. 

Multigrid solver based on DGS smoother converges uniformly. 

To-Do: test other boundary conditions.
