function showsolutionP2(node,elem,u)
%% SHOWSOLUTIONCR plots the CR function u on a triangular mesh in 2-D.

[node,elem] = uniformrefine(node,elem);
showsolution(node,elem,u);