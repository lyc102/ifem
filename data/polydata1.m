function pde = polydata1
%% POLYDATA1 polynomial data for Poisson equation
%
% f = -6*(x+y);
% u = x^3 + y^3;
% Du = (3x^2, 3y^2);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du);

    % load data (right hand side function)
    function z = f(p)
    x = p(:,1); y = p(:,2);
    z = -6*(x+y);
    end
    % Dirichlet boundary condition
    function z = g_D(p)
    z = exactu(p);
    end
    % exact solution
    function z = exactu(p)
    x = p(:,1); y = p(:,2);
    z = x.^3 + y.^3;
    end
    % Derivative of the exact solution
    function z = Du(p)
    x = p(:,1); y = p(:,2);
    z(:,1) = 3*x.^2;
    z(:,2) = 3*y.^2;
    end
end