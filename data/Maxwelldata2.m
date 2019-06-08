function pde = Maxwelldata2
%% MAXWELLDATA2 non-homogenous Dirichlet/Neumann boundary condition
%
%       curl(mu^(-1)curl u) -     epsilon u = J    in \Omega,  
%                                n \corss u = n � g_D  on \Gamma_D,
%                  n \cross (mu^(-1)curl u) = n � g_N  on \Gamma_N.
%
%   mu = 1; 
%   epsilon = - 1; % positive type
%   u = [0 cos(x) cos(x)];
%   curlu = [0 sin(x) -sin(x)];
%   J = 2*u;
%   g_N = n � curl u;
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
        s = 2*exactu(p);
    end

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