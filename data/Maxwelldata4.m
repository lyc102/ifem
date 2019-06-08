function pde = Maxwelldata4
%% MAXWELLDATA4 homogenous Dirichlet boundary condition
%
%       curl(mu^(-1)curl u) - epsilon u = J    in \Omega,  
%                                 n × u = n × g_D  on \Gamma_D,
%                   n × (mu^(-1)curl u) = n × g_N  on \Gamma_N?
%
%   mu = 1; 
%   epsilon = - 1; % positive type
%   u = [(x^2-1)*(y^2-1)*(z^2-1) 0 0];
%   curlu = [0, 2*(x^2-1)*(y^2-1)*z, -2*(x^2-1)*(z^2-1)*y];;
%   curlcurl u = [-2*(x^2-1)*(z^2-1)-2*(x^2-1)*(y^2-1), ...
%                   4*(z^2-1)*x*y, 4*x*z*(y^2-1)];
%   J = curlcul u + u;
%
%  - u satisfies zero Dirichlet boundary condition at cube [-1,1]^3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.mu = 1;
pde.epsilon = -1;
pde.J = @J;
pde.exactu = @exactu;
pde.g_D = @exactu;
pde.curlu = @curlu;
pde.g_N = @g_N;

    function s = J(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s = [-2*(x.^2-1).*(z.^2-1)-2*(x.^2-1).*(y.^2-1), ...
         4*(z.^2-1).*x.*y, 4*x.*z.*(y.^2-1)];
        s = s + exactu(p);
    end

    function s = exactu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s = [(x.^2-1).*(y.^2-1).*(z.^2-1), 0*y, 0*z];
    end

    function s = curlu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s = [0*x, 2*(x.^2-1).*(y.^2-1).*z, -2*(x.^2-1).*(z.^2-1).*y];
    end

    function s = g_N(p,n)
        s = curlu(p);
        s = cross(n,s,2);
    end
end