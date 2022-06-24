%% 
%
% Created by Jie Zhou and clean up by Long Chen.
%
% Reference
%
% https://fenicsproject.org/olddocs/dolfin/1.3.0/python/demo/documented/cahn-hilliard/python/documentation.html
%  https://fenicsproject.org/docs/dolfin/1.3.0/python/demo/documented/cahn-hilliard/python/documentation.html

close all
%% Parameters
%  maxN = 2e3;     theta = 0.5;    maxIt = 1; 
%  N = zeros(maxIt,1);     errL2 = zeros(maxIt,1);     errH1 = zeros(maxIt,1);

%%  Generate an initial mesh
h = 1/80;
[node,elem] = squaremesh([0,1,0,1],h);

%% Parametere
pde = CahnHilliarddata;
global efsilonsquare 
efsilonsquare = 1.0e-2;
tao  = 5.0e-6;  %time step
times = 101;     

%% initial condition;
u0 = pde.initial_u(node); % initial vavlue u
w0 = pde.initial_w(node); % initial vavlue w

%% Solve
[w0,u0] = CahnHilliardP1(node,elem,pde,w0,u0,tao,times);