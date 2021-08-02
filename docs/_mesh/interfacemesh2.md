# Interface-fitted Mesh Generator: Two Dimensions

 

We explain a simple mesh generator for generating an interface-fitted mesh in two dimensions. 

Consider a square domain $\Omega$ with an interface $\Gamma$ in it:

<img src="../assets/images/mesh/interfacedomain.png" alt="interfacedomain" style="zoom:67%;" />



**Algorithm:  2D Interface-fitted Mesh Generation**

Grid size: h; Level set function, φ(x); Square domain, Ω; OUTOUT:

Interface-fitted mesh T ;

1. 1:  Find the cut points, the Cartesian mesh points near or on

   the interface and some auxiliary points.

2. 2:  Construct a Delaunay triangulation on these points.

3. 3:  Post processing. Remove triangles away from the interface

   and merge all uncut elements