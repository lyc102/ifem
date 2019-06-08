function pde = fracLapdata2
%% FRACLAPDATA2 f = sin(m*pi*x1).*sin(n*pi*x2)
%
% s a parameter in (0,1)
% m,n given integers
% lambda_(m,n) = pi^2*(m^2+n^2)
% sqrtlambda = sqrt((m^2+n^2))*pi
% u = 2^(1-s)/gamma(s)*(sqrtlambda*y)^s*besselk(s,sqrtlambda*y).*sin(m*pi*x1)*sin(n*pi*x2);
% u|_(y==0) =  sin(m*pi*x1)*sin(n*pi*x2);
% f = (-Delta)^s u = lambda^s u 
% g_N = ds*f;
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'Du',@Du);

    % exact solution
    function u =  exactu(p) 
        global s
        m = 2; n = 2;
        u = zeros(size(p,1),1);
        x1 = p(:,1); x2 = p(:,2); 
        if size(p,2) == 3
            y = p(:,3);
        else
            y = zeros(size(p,1),1);
        end
        idx = (y>eps);
        sqrtlambda = sqrt(m^2+n^2)*pi;
        z = sqrtlambda*y(idx);
        C = 2^(1-s)/gamma(s);
        lambda = pi^2*(m^2+n^2);
        u(idx) =  lambda^(-s)*C*z.^s.*besselk(s,z).*sin(m*pi*x1(idx)).*sin(n*pi*x2(idx));
        u(~idx) = lambda^(-s)*sin(m*pi*x1(~idx)).*sin(n*pi*x2(~idx));  % at y = 0.  %% NOS 2013 p.36
    end
    function z = f(p)
%         global s
        m = 2; n = 2;
        x1 = p(:,1); x2 = p(:,2); % y = p(:,3);
        z = sin(m*pi*x1).*sin(n*pi*x2);
    end
    function uprime =  Du(p)
        % not right. modify later.
        global s
        m = 2; n = 2;
        x1 = p(:,1); x2 = p(:,2); y = p(:,3);
        sqrtlambda = sqrt(m^2+n^2)*pi;
        z = sqrtlambda*y;
        C = 2^(1-s)/gamma(s);
        uprime(:,1) = C*m*pi*z.^s.*besselk(s,z).*cos(m*pi*x1).*sin(n*pi*x2);
        uprime(:,2) = C*n*pi*z.^s.*besselk(s,z).*sin(m*pi*x1).*cos(n*pi*x2);
        uprime(:,3) = -C*sqrtlambda*z.^s.*besselk(1-s,z).*sin(m*pi*x1).*sin(n*pi*x2);
    end
end