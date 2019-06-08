function pde = sincosdata3test
%% SINCOSDATA3 trigonometric data for Poisson equation in 3-D
%
%     f = 3*pi^2*cos(pi*x)*cos(pi*y)*cos(pi*z);
%     u = cos(pi*x)*cos(pi*y)*cos(pi*z);
%     Du = (-pi*sin(pi*x)*cos(pi*y)*cos(pi*z), 
%           -pi*cos(pi*x)*sin(pi*y)*cos(pi*z),
%           -pi*cos(pi*x)*cos(pi*y)*sin(pi*z));
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du);

    % load data (right hand side function)
    function s = f(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = 3*pi^2*cos(pi*x).*cos(pi*y).*cos(pi*z);
    end
    % exact solution
    function s = exactu(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = cos(pi*x).*cos(pi*y).*cos(pi*z); 
    end
    % Dirichlet boundary condition
    function s = g_D(p)
    s = exactu(p);
    end

    % Derivative of the exact solution
    function s = Du(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s(:,1) = -pi*sin(pi*x).*cos(pi*y).*cos(pi*z);
    s(:,2) = -pi*cos(pi*x).*sin(pi*y).*cos(pi*z);
    s(:,3) = -pi*cos(pi*x).*cos(pi*y).*sin(pi*z);
    end
end