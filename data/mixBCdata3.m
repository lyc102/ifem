function pde = mixBCdata3
%% MIXBCDATA3 mix boundary condition data for Poisson equation in 3-D
%
%     f = 3*sin(x)*sin(y)*sin(z);
%     u = sin(x)*sin(y)*sin(z);
%     Du = (cos(x)*sin(y)*sin(z),
%           sin(x)*cos(y)*sin(z),
%           sin(x)*sin(y)*cos(z));
%
% Created by Ming Wang.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'g_N',@g_N,'Du',@Du);

    % load data (right hand side function)
    function s = f(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = 3*sin(x).*sin(y).*sin(z);
    end
    % exact solution
    function s = exactu(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = sin(x).*sin(y).*sin(z);
    end
    % Dirichlet boundary condition
    function s = g_D(p)
    s = exactu(p);
    end
    % Neumann boundary condition,[-1,1]*[-1,1]*[-1,1].
    eps = 1.0e-14; % a small quantity.
    function s = g_N(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    x0 = abs(x+1)< eps; x1 = abs(x-1)< eps;
    y0 = abs(y+1)< eps; y1 = abs(y-1)< eps;
    z0 = abs(z+1)< eps; z1 = abs(z-1)< eps;
    % normal derivative at corners doesn't make sense.
    s = x0.*(-cos(x).*sin(y).*sin(z))+x1.*(cos(x).*sin(y).*sin(z))+ ...
        y0.*(-sin(x).*cos(y).*sin(z))+y1.*(sin(x).*cos(y).*sin(z))+ ...
        z0.*(-sin(x).*sin(y).*cos(z))+z1.*(sin(x).*sin(y).*cos(z));
    end
    % Derivative of the exact solution
    function s = Du(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s(:,1) = cos(x).*sin(y).*sin(z);
    s(:,2) = sin(x).*cos(y).*sin(z);
    s(:,3) = sin(x).*sin(y).*cos(z);
    end
end