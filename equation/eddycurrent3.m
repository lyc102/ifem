function [u,edge,eqn,info] = eddycurrent3(node,elem,bdFlag,pde,option,varargin)
%% EDDYCURRENT3 Eddy current equation in 3-D using lowest order edge element.
%
% [u,edge,eqn,info] = EDDYCURRENT3(node,elem,bdFlag,pde,option,HB)
% solves
%
%   curl(mu^(-1)curl u) + beta*u = J    in Omega,
%                       n x u = n x g_D on Gamma_D,
%         n x (mu^(-1)curl u) = n x g_N on Gamma_N.
%
%
% Input fields in pde:
%   - pde.mu, pde.beta, pde.J, pde.g_D, pde.g_N
% Example
%
% See also eddycurrent
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
if ~exist('pde','var') || isempty(pde), pde = struct(); end

%% Defaults for PDE data
if ~isfield(pde,'mu'), pde.mu = 1; end
if ~isfield(pde,'beta'), pde.beta = 1; end
if ~isfield(pde,'J'), pde.J = []; end
if ~isfield(pde,'g_D'), pde.g_D = []; end
if ~isfield(pde,'g_N'), pde.g_N = []; end

%% Map eddy current coefficients to Maxwell coefficients
pdeMaxwell = pde;
pdeMaxwell.omega = 1;
if isnumeric(pde.beta)
    pdeMaxwell.epsilon = -pde.beta;
else
    pdeMaxwell.epsilon = @(x) -pde.beta(x);
end

%% Solver option compatibility
if isfield(option,'solver')
    if strcmp(option.solver,'nosolve') || strcmp(option.solver,'notsolve')
        option.solver = 'none';
    end
end

%% Delegate to 3-D Maxwell solver
if nargin >= 6
    HB = varargin{1};
    [u,edge,eqn,info] = Maxwell(node,elem,bdFlag,pdeMaxwell,option,HB);
else
    [u,edge,eqn,info] = Maxwell(node,elem,bdFlag,pdeMaxwell,option);
end

%% Convert returned mass term to eddy-current sign convention
if ~isempty(eqn) && isfield(eqn,'M')
    eqn.M = -eqn.M;
end

end
