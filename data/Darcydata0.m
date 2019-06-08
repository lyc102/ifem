function pde = Darcydata0
%% Darcydata0 trigonometric  data for Poisson equation 
%
%     f = 8*pi^2*sin(2*pi*x)*cos(2*pi*y); f = -\Delta p
%     p = sin(2*pi*x)*cos(2*pi*y);
%     Dp = (2*pi*cos(2*pi*x)*cos(2*pi*y), -2*pi*sin(2*pi*x)*sin(2*pi*y));
%     u = Dp;
%     du/dn = g_N on [0,1]^2.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactp',@exactp,'exactu',@exactu,'g_D',@g_D,'g_N',@g_N);

    % load data (right hand side function)
    function rhs =  f(p)
        x = p(:,1); y = p(:,2);
        rhs =  8*pi^2*sin(2*pi*x).*cos(2*pi*y);
    end
    % exact solution
    function s =  exactp(pt)
        x = pt(:,1); y = pt(:,2);
        s = sin(2*pi*x).*cos(2*pi*y);
    end
    % the derivative of the exact solution
    function s = exactu(p)
        x = p(:,1); y = p(:,2);
        s(:,1) = 2*pi*cos(2*pi*x).*cos(2*pi*y);
        s(:,2) = -2*pi*sin(2*pi*x).*sin(2*pi*y);
    end
    function u = g_D(p)
        u = exactp(p);
    end
    % Neumann boundary condition 
    function f = g_N(p,vargin)
        if nargin > 1
            f = dot(exactu(p),vargin,2);
        else
            f = zeros(size(p,1),1);
            x = p(:,1); y = p(:,2);
            uprime = exactu(p);
            leftbd = (abs(x)<eps);  % n = (-1,0); 
            f(leftbd) = - uprime(leftbd,1);
            rightbd = (abs(x-1)<eps); % n = (1,0); 
            f(rightbd) = uprime(rightbd,1);
            topbd = (abs(y-1)<eps);   % n = (0,1)
            f(topbd) = uprime(topbd,2);
            bottombd = (abs(y)<eps);% n = (0,-1)
            f(bottombd) = - uprime(bottombd,2);    
        end
    end
end