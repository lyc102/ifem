function pde = Helmholtzdata1
%% HELMHOTLZDATA1 trigonometric  data for Helmholtz equation
%
%     f = (5*pi^2 - k^2)*sin(pi*x)*sin(2*pi*y);
%     u = sin(pi*x)*sin(2*pi*y);
%   g_D = 0
%
% Reference: Example 6.1 Closed-off problem in 
% Erlangga, Y. A., Vuik, C., & Oosterlee, C. W. (2004). On a
% class of preconditioners for solving the Helmholtz equation. Applied
% Numerical Mathematics, 50(3-4), 409?425.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'k2',@k2,'Du',@Du);

    % load data (right hand side function)
    function rhs =  f(p)
        global k
        x = p(:,1); y = p(:,2);
        rhs = (-k^2 + 5*pi^2)*sin(pi*x).*sin(2*pi*y);
    end
    % exact solution
    function u =  exactu(p)
        x = p(:,1); y = p(:,2);
        u = sin(pi*x).*sin(2*pi*y);
    end
    % Derivative of the exact solution
    function uprime =  Du(p)
        x = p(:,1); y = p(:,2);
        uprime(:,1) = pi*cos(pi*x).*sin(2*pi*y);
        uprime(:,2) = 2*pi*sin(pi*x).*cos(2*pi*y);
    end
    % wave number
    function wavenumber = k2(~)
        global k
        wavenumber = k.^2;
    end
end