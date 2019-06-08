function pde = Maxwelldata3
%% MAXWELLDATA3 polynomial data
%       curl(mu^(-1)curl u) - epsilon u = J    in \Omega,  
%                                 n × u = n × g_D  on \Gamma_D,
%                   n × (mu^(-1)curl u) = n × g_N  on \Gamma_N?
%
%   mu = 1; 
%   epsilon = - 1; % positive type
%
%     u = grad (x^n+y^n+z^n) = [n*x.^(n-1), n*y.^(n-1), n*z.^(n-1)];
%     curlu = 0;
%     J = u;
%
%  u satisfies the Neumann boundary condition curlu cross n = 0.
%  curl u = 0 equivalent to L2 projection.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.mu = 1;
pde.epsilon = -1;
pde.J = @J;
pde.exactu = @exactu;
pde.g_D = @exactu;
pde.curlu = [0 0 0];

    function s = J(p)
    s = exactu(p);
    end

    function s = exactu(p)
    n = 4;
    x = p(:,1); y = p(:,2); z = p(:,3);
    % u = grad (x^n+y^n+z^n);
    s = [n*x.^(n-1), n*y.^(n-1), n*z.^(n-1)];
    % curl u = 0;
    end
end