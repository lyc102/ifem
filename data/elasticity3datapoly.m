function pde = elasticity3datapoly(para)
%% ELASTICITYDATA3 data for the elasticity problem in three dimensions
%
% Modified from elasticitydata by Huayi Wei.
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
    function s = f(p)
% f = - mu \Delta u - (mu + lambda) grad(div u) 
       x = p(:,1); y = p(:,2); z = p(:,2);
       f1 = lambda*(32*y.*z.*(y - 1).*(z - 1) + 64*x.*y.*(2*z - 1).*(y - 1) ...
           + 32*x.*z.*(2*y - 1).*(z - 1) + 64*y.*(2*z - 1).*(x - 1).*(y - 1) ...
           + 32*z.*(2*y - 1).*(x - 1).*(z - 1)) + 32*mu*z.*(z - 1).*(4*x.*y - 2*y - 3*x + x.^2 + 1) ...
           + 32*mu*y.*(y - 1).*(8*x.*z - 4*z - 5*x + x.^2 + 2) + 64*mu*y.*z.*(y - 1).*(z - 1);
       f2 = lambda*(64*x.*z.*(x - 1).*(z - 1) + 64*x.*y.*(2*z - 1).*(x - 1) ...
           + 16*y.*z.*(2*x - 1).*(z - 1) + 64*x.*(2*z - 1).*(x - 1).*(y - 1) ...
           + 16*z.*(2*x - 1).*(y - 1).*(z - 1)) + 64*mu*x.*(x - 1).*(4*y.*z - 2*z - 3*y + y.^2 + 1) ...
           + 16*mu*z.*(z - 1).*(4*x.*y - 6*y - 2*x + 4*y.^2 + 1) + 128*mu*x.*z.*(x - 1).*(z - 1);
       f3 = lambda*(128*x.*y.*(x - 1).*(y - 1) + 32*x.*z.*(2*y - 1).*(x - 1) ...
           + 16*y.*z.*(2*x - 1).*(y - 1) + 32*x.*(2*y - 1).*(x - 1).*(z - 1) ...
           + 16*y.*(2*x - 1).*(y - 1).*(z - 1)) + 32*mu*x.*(x - 1).*(4*y.*z - 6*z - 2*y + 4*z.^2 + 1) ...
           + 16*mu*y.*(y - 1).*(4*x.*z - 10*z - 2*x + 8*z.^2 + 1) + 256*mu*x.*y.*(x - 1).*(y - 1);
       s = [f1 f2 f3];            
    end
    function u = exactu(p)
       x = p(:,1); y = p(:,2); z = p(:,2);
       u1 = 16*x.*(1-x).*y.*(1-y).*z.*(1-z);
       u2 = 32*x.*(1-x).*y.*(1-y).*z.*(1-z);
       u3 = 64*x.*(1-x).*y.*(1-y).*z.*(1-z);
       u = [u1 u2 u3];
    end
    function s = g_D(p)
       s = exactu(p); 
    end
end