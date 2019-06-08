function pde = fracLapdata1
%% FRACLAPDATA1 data for fractional Laplacian problem
%
% s a parameter in (0,1)
% k is given integer
% u = 2^(1-s)/gamma(s)*(k*pi*y)^s*besselk(s,k*pi*y).*sin(k*pi*x);
% u|_(y==0) =  sin(k*pi*x);
% g_N = ds*(k*pi)^(2*s)*sin(k*pi*x);
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',0,'exactu',@exactu,'g_D',0,'g_N',@g_N,'d',@d,'Du',@Du);

    % exact solution
    function u =  exactu(p)
        global s
        u = zeros(size(p,1),1);
        k = 3;
        x = p(:,1); y = p(:,2);
        idx = (y>eps);
        sqrtlambda = k*pi;
        z = sqrtlambda*y(idx);
        C = 2^(1-s)/gamma(s);
        u(idx) =  C*z.^s.*besselk(s,z).*sin(sqrtlambda*x(idx));
        u(~idx) = sin(sqrtlambda*x(~idx));  % at y = 0.  %% NOS 2013 p.31
    end
    function z = g_N(p)
        global s
        k = 3;
        x = p(:,1); %y = p(:,2);
        ds = 2^(1-2*s)*gamma(1-s)/gamma(s);  % (2.23) in NOS 2013
        sqrtlambda = k*pi;
        z = ds*sqrtlambda^(2*s)*sin(sqrtlambda*x);
    end
    function uprime =  Du(p)
        global s
        k = 3;
        x = p(:,1); y = p(:,2);
        sqrtlambda = k*pi;
        z = sqrtlambda*y;
        uprime(:,1) = 2^(1-s)/gamma(s)*sqrtlambda*z.^s.*besselk(s,z).*cos(sqrtlambda*x);
        uprime(:,2) = -2^(1-s)/gamma(s)*sqrtlambda*z.^s.*besselk(1-s,z).*sin(sqrtlambda*x);
    end
end