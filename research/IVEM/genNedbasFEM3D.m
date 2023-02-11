function fem = genNedbasFEM3D(mesh,type)

%% Usage: Form Conforming P1 FEM degrees of freedom structure
%
% INPUTS:
% mesh --- a struct data contains very rich mesh information.
%
% OUTPUTS:
% fem --- a struct data contains the following fields:
%         fem.p: (x,y,z) coordinate of each vertex w.r.t a global DoF
%         fem.t: indices of global DoF in each element
%         fem.type: type of FEM
%         fem.ldof: number of local DoF on each element
%         fem.bas: three dimension matrix structure contains
%                  basis function on each element

% Last Modified: 07/02/2020 by Xu Zhang

%% 1. Form p,t
p = mesh.p; t = mesh.t; 

%% 2. Form basis in vector form
bas = bas3DP1(p,t);
if type == 1
    bas = bas3D_ned1(bas);
end
fem = struct('p',p, 't',t, 'type','P1', 'ldof',6, 'bas',bas); 