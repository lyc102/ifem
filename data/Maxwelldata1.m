function pde = Maxwelldata1
%% MAXWELLDATA1 positive definite and homogenous Neumann boundary condition
%
%  curl(mu^(-1)curl u) - omega^2 epsilon u = J    in \Omega,  
%                                    n × u = n × g_D  on \Gamma_D,
%                      n × (mu^(-1)curl u) = n × g_N  on \Gamma_N?
%
%  mu = 1; 
%  omega = 1;
%  epsilon = - 1; positive type
%  u = grad(sin(x) + cos(y) + z^2+z) = [cos(x), -sin(y), 2*z+1];
%  curlu = 0;
%  J = u;
%     
%  * curl u = 0 equivalent to L2 projection.
%  * u satisfies the zero Neumann boundary condition curlu cross n = 0.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.mu = 1;
pde.epsilon = -1;
pde.J = @J;
pde.exactu = @exactu;
pde.g_D = @exactu;
pde.curlu = [0 0 0];
pde.g_N = 0;

    function s = J(p)
    s = exactu(p);
    end

    function s = exactu(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = [cos(x), -sin(y), 2*z+1];
    end
end