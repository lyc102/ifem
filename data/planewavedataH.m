function pde = planewavedataH
%% PLANEWAVEDATAH plane wave solution of Maxwell equations in terms of H
%
%   epsilon = mu = omega = 1;
%   E = exp(-i*omega*x*d)*P;
%   H = -exp(-i*omega*x*d)*d\times P
%   where d\in S^2, P \in C^3 and d*P = 0.
%   J = 0;
%   g_N = curl H /epsilon = -i*omega*E.
%
% For Maxwell equation in terms of E, we set u = E.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.mu = 1;
pde.epsilon = 1;
pde.J = 0;
pde.exactu = @H;
pde.g_D = @H;
pde.curlu = @E;
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
        s = -li*omega*E(x);
        s = cross(n,s,2);
    end
end