function pde = planewavedata1
%% PLANEWAVEDATA1 plane wave solution of Maxwell equations: real solution.
%
%       curl(curl u) - u = 0    in \Omega,  
%                  n × u = n × g_D  on \Gamma_D,
%             n × curl u = n × g_N  on \Gamma_N?
%
%   mu = epsilon = omega = 1;
%   u = [0 cos(x) cos(x)];
%   curlu = [0 sin(x) -sin(x)];
%   curl curl u = u;
%   J = 0;
%   g_N = curl u;
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.mu = 1;
pde.epsilon = 1;
pde.J = 0;
pde.exactu = @exactu;
pde.g_D = @exactu;
pde.curlu = @curlu;
pde.g_N = @g_N;

    function s = exactu(p)
        x = p(:,1); %y = p(:,2); z = p(:,3);
        s = [0*x, cos(x), cos(x)];
    end

    function s = curlu(p)
        x = p(:,1); %y = p(:,2); z = p(:,3);
        s = [0*x, sin(x), -sin(x)];
    end

    function s = g_N(p,n)
        s = curlu(p);
        s = cross(n,s,2);
    end
end