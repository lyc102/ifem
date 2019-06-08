function pde = planewavedata
%% PLANEWAVEDATA plane wave solution of Maxwell equations: complex solution
%
%   epsilon = mu = 1;
%   E = exp(-i*omega*x*d)*P;
%   H = -exp(-i*omega*x*d)*d\times P
%   where d\in S^2, xi \in C^3 and d\times xi = 0.
%   J = 0;
%
% For Maxwell equation in terms of E, we set u = E.
%       curl(curl u) - u = 0    in \Omega,  
%                  n × u = n × g_D  on \Gamma_D,
%             n × curl u = n × g_N  on \Gamma_N?
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.mu = 1;
pde.epsilon = 1;
pde.J = 0;
pde.exactu = @E;
pde.g_D = @E;
pde.curlu = @H;
pde.g_N = @g_N;

li = sqrt(-1);
    function s = H(x)
        global d P omega
        s = -exp(-1i*omega*(x(:,1)*d(1)+x(:,2)*d(2)+x(:,3)*d(3)))*cross(d,P);
    end

    function s = E(x)
        global d P omega
        s = exp(-1i*omega*(x(:,1)*d(1)+x(:,2)*d(2)+x(:,3)*d(3)))*P;
    end

    function s = g_N(x,n)
        global omega
        s = li*omega*H(x);
        s = cross(n,s,2);
    end
end