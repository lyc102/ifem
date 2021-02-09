function pde = quadCurlDataSmooth1
% (curl)^4 u = f on \Omega = [0,pi]^3
% u x n = 0
% (curl u) x n = 0
% div u = g
%
pde = struct('exactu', @exactu,'g_D', @g_D, 'g', @divu, ...
             'curlu', @curlu, 'curlcurlu',@curlcurlu, ...
             'tricurlu',@tricurlu,'quadcurlu',@quadcurlu, 'pentacurlu',@pentacurlu);

Scaling = 10; % scaling
         
    function s = exactu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s(:,1) = 0*x; 
        s(:,2) = 0*y;
        s(:,3) = (sin(x)).^2.*(sin(y)).^2.*sin(z);
        s = s/Scaling;
    end

    function s = curlu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s(:,1) =  sin(2*y).*(sin(x)).^2.*sin(z);
        s(:,2) = -sin(2*x).*(sin(y)).^2.*sin(z);
        s(:,3) = 0*z;
        s = s/Scaling;
    end

    function s = curlcurlu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s(:,1) = sin(2*x).*(sin(y)).^2.*cos(z);
        s(:,2) = (sin(x)).^2.*sin(2*y).*cos(z);
        s(:,3) = 4*(sin(x)).^2.*(sin(y)).^2.*sin(z) ...
               - 2*(cos(y)).^2.*(sin(x)).^2.*sin(z) ...
               - 2*(cos(x)).^2.*(sin(y)).^2.*sin(z);
        s = s/Scaling;
    end


    function s = tricurlu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s(:,1) = 7*(sin(x)).^2.*sin(2*y).*sin(z) - 2*(cos(x)).^2.*sin(2*y).*sin(z);
        s(:,2) = 2*(cos(y)).^2.*sin(2*x).*sin(z) - 7*sin(2*x).*(sin(y)).^2.*sin(z);
        s(:,3) = 0*z;
        s = s/Scaling;
    end
    

    function s = quadcurlu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s(:,1) = 7*cos(z).*sin(2*x).*(sin(y)).^2 - 2*(cos(y)).^2.*cos(z).*sin(2*x);
        s(:,2) = 7*cos(z).*(sin(x)).^2.*sin(2*y) - 2*(cos(x)).^2.*cos(z).*sin(2*y);
        s(:,3) = 8*(cos(x)).^2.*(cos(y)).^2.*sin(z) - 18*(cos(x)).^2.*(sin(y)).^2.*sin(z) ...
              - 18*(cos(y)).^2.*(sin(x)).^2.*sin(z) + 28*(sin(x)).^2.*(sin(y)).^2.*sin(z);
        s = s/Scaling;
    end

    function s = pentacurlu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s(:,1) = 53*(sin(x)).^2.*sin(2*y).*sin(z) - 28*(cos(x)).^2.*sin(2*y).*sin(z);
        s(:,2) = 28*(cos(y)).^2.*sin(2*x).*sin(z) - 53*sin(2*x).*(sin(y)).^2.*sin(z);
        s(:,3) = 0*z;
        s = s/Scaling;
    end

    function s = g_D(p) % Dirichlet boundary condition for u
        s = exactu(p);
    end

    function s = divu(p) 
        x = p(:,1); y = p(:,2); z = p(:,3);
        s = (sin(x)).^2.*(sin(y)).^2.*cos(z);
        s = s/Scaling;
    end
end