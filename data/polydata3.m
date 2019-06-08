function pde = polydata3
%% POLYDATA3 polynomial data for Poisson equation in 3-D
%
% f = -6*(x + y + z);
% u = x^3 + y^3 + z^3;
% Du = (3x^2, 3y^2, 3z^2);

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du);

    % load (right hand side function)
    function s = f(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = -6*(x + y + z);
    end
    % Dirichlet boundary condition
    function s = g_D(p)
    s = exactu(p);
    end
    % exact solution
    function s = exactu(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = x.^3 + y.^3 + z.^3;
    end
    % Derivative of the exact solution
    function s = Du(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s(:,1) = 3*x.^2;
    s(:,2) = 3*y.^2;
    s(:,3) = 3*z.^2;
    end
end