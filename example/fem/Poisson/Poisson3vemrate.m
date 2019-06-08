clear all; close all;
pde = sincosdata3;
option.h0 = 0.1;
option.maxIt = 3;
option.solver = 'amg';
cube = [-1,1,-1,1,-1,1];
vemPoisson3(cube, pde, option);