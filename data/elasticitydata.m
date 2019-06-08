function pde = elasticitydata(para)
%% ELASTICITYDATA data for elasticity problem
% 
% - mu \Delta u - (mu + lambda) grad(div u) = f    in \Omega
%                                         u = g_D  on \partial \Omega
%
% Created by Huayi Wei Monday, 27 June 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Lame constants
if nargin == 0
    lambda = 1;
    mu = 1;
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'lambda') || isempty(para.lambda)
        lambda = 1;
    else
        lambda = para.lambda;
    end
    if ~isfield(para,'mu') || isempty(para.mu)
        mu = 1;
    else
        mu = para.mu;
    end
end

pde = struct('lambda',lambda,'mu',mu, 'f', @f, 'exactu',@exactu,'g_D',@g_D);
%%%%%% subfunctions %%%%%%
    function z = f(p)
       z = mu*2*pi^2.*exactu(p);
    end
    function z = exactu(p)
       x = p(:,1); y = p(:,2);
       z = [cos(pi*x).*cos(pi*y), sin(pi*x).*sin(pi*y)];    
    end
    function z = g_D(p)
       z = exactu(p); 
    end
end