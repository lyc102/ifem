function [theta,coarsen,interpolation,preconditioner,mu] = amgoption(option)
%% AMGOPTIONS default options for algebraic multigrid solver
%  
%   theta = 0.025;	cmethod = 'C'; interpolation = 'S'; preconditioner = 'W';
%  
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.


if isfield(option,'theta')
    theta = option.theta;
else
    theta = 0.025;
end

if isfield(option,'coarsen')
    coarsen = option.coarsen;
else
    coarsen = 'C';
end

if isfield(option,'interpolation')
    interpolation = option.interpolation;
else
    interpolation = 'T';
end
if isfield(option,'preconditioner')
    preconditioner = option.preconditioner;
else
    preconditioner = 'W';
end
if isfield(option,'smoothingstep')  % smoothing steps
    mu = option.smoothingstep;
else
    mu = 2;
end