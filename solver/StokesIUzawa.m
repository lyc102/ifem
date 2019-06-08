function [u,p] = StokesIUzawa(u,p,f,g,A,B,auxMat,elem,smootherOpt)
%% STOKESIUzawa Inexact Uzawa smoother for the Stokes eqns.
%
%  In matrix form, we need to solve the equations
%
%         Lx =  |A B'| |u| = |f|
%               |B 0 | |p| = |0|
%
%  Define
%            M = |I -S_u^-1 B'|   T = |S_u    0|
%                |     I      |       |B   -S_p|
%  The matrix form DGS update can be written as
%
%    |uk+1|   |uk|            |ru|     with ru = f-Auk -B'pk,
%    |    | = |  | + M*inv(T)*|  |
%    |pk+1|   |pk|            |rp|     and  rp = 0-Buk
%
%  Created by Ming Wang (with discussion of Long Chen) at Jan, 2012. 
%  Revised at Sept 2012.
%

%% Parameters
itStep = smootherOpt.itStep;
Su = auxMat.DA; % DA = 2*diag(A)
Bt = auxMat.Bt;
BinvDABt = auxMat.BinvDABt; % invDA = 1/(2*diag(A));
% mg options
mgoption.printlevel = 0; mgoption.maxIt = 1; mgoption.N0=50;
mgoption.solver = 'Vcycle'; 
mgoption.mu = 2;

%% DGS relaxation step
for k = 1: itStep
    % Step 1: relax Momentum eqns
    if isempty(u) && isempty(p)
        ru = f;
        u = zeros(size(f)); p = zeros(size(g));
    else
        ru = f-Bt*p-A*u;
    end
    u = u + ru./Su;
    rp = g-B*u;
    rp = rp - mean(rp);
    % Step 2: relax transformed Continuity eqns
    dp = mg(BinvDABt,-rp,elem,mgoption);
    % Step 3: transform the correction back to the original variables.
    p = p + dp;
    u = u - (Bt*dp)./Su;
end
end