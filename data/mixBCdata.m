function pde = mixBCdata
%% MIXBCDATA mix boundary condition data for Poisson equation
%     f = 2*sin(x)*sin(y);
%     u = sin(x)*sin(y);
%     Du = (cos(x)*sin(y), sin(x)*cos(y));
%
%
% Created by Ming Wang.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'g_N',@g_N,'Du',@Du,'d',@d);
pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'g_N',@g_N,'Du',@Du,'Du1',@Du1,'Du2',@Du2);
    % load data (right hand side function)
    function z = f(p)
    x = p(:,1); y = p(:,2);
    z = 2*sin(x).*sin(y);
    end
    % exact solution
    function z = exactu(p)
    x = p(:,1); y = p(:,2);
    z = sin(x).*sin(y);
    end
    % Dirichlet boundary condition
    function z = g_D(p)
    z = exactu(p);
    end
    % Neumann boundary condition,[0,1]*[0,1]
    function z = g_N(p)
    x = p(:,1); y = p(:,2);
    x0 = abs(x)< eps; x1 = abs(x-1)< eps;
    y0 = abs(y)< eps; y1 = abs(y-1)< eps;
    % normal derivative at corners doesn't make sense.
    z = x0.*(-cos(x).*sin(y))+x1.*(cos(x).*sin(y))+ ...
        y0.*(-sin(x).*cos(y))+y1.*(sin(x).*cos(y));
    end
    % Derivative of the exact solution
    function z = Du(p)
    x = p(:,1); y = p(:,2);
    z(:,1) = cos(x).*sin(y);
    z(:,2) = sin(x).*cos(y);
    end

    function z = Du1(p) 
    x = p(:,1); y = p(:,2);
    z = cos(x).*sin(y);
    end

    function z = Du2(p)
    x = p(:,1); y = p(:,2);
    z = sin(x).*cos(y);
    end
    % function z = d(p)
    % z = 1/2;
    % end
end