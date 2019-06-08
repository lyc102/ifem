function pde = sincosdata
%% SINCOSDATA trigonometric  data for Poisson equation
%
%     f = 2*pi^2*cos(pi*x)*cos(pi*y);
%     u = cos(pi*x)*cos(pi*y);
%     Du = (-pi*sin(pi*x)*cos(pi*y), -pi*cos(pi*x)*sin(pi*y));
%
% The u satisfies the zero flux condition du/dn = 0 on boundary of [0,1]^2
% and thus g_N is not assigned.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du);

    % load data (right hand side function)
    function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs =  2*pi^2*cos(pi*x).*cos(pi*y);
    end
    % exact solution
    function u =  exactu(p)
    x = p(:,1); y = p(:,2);
    u =  cos(pi*x).*cos(pi*y);
    end
    % Dirichlet boundary condition
    function u =  g_D(p)
    u =  exactu(p);
    end
    % Derivative of the exact solution
    function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1) = -pi*sin(pi*x).*cos(pi*y);
    uprime(:,2) = -pi*cos(pi*x).*sin(pi*y);
    end
end