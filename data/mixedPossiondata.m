function pde = mixedPossiondata
%% MIXBCDATA3 mix boundary condition data for Poisson equation in 3-D
%
%     u = sin(pi*x)^2*sin(pi*y)^2*sin(pi*z)^2;
%     f = - div Du
%
% Created by Yongke Wu.

pde.Du = @Du;
pde.f = @f;
pde.exactu = @u;

    function s = u(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s = (sin(pi.*x)).^2.*(sin(pi.*y)).^2.*(sin(pi.*z)).^2;
    end

    function s = Du(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s = 2*pi.*[sin(pi.*x).*cos(pi.*x).*(sin(pi.*y)).^2.*(sin(pi.*z)).^2,...
                   (sin(pi.*x)).^2.*sin(pi.*y).*cos(pi.*y).*(sin(pi.*z)).^2,...
                   (sin(pi.*x)).^2.*(sin(pi.*y)).^2.*sin(pi.*z).*cos(pi.*z)];
    end

    function s = f(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        s = - 2*pi^2*((cos(pi*x)).^2.*(sin(pi*y)).^2.*(sin(pi*z)).^2 ...
                  + (cos(pi*y)).^2.*(sin(pi*x)).^2.*(sin(pi*z)).^2 ...
                  + (cos(pi*z)).^2.*(sin(pi*x)).^2.*(sin(pi*y)).^2 ...
                  - 3*(sin(pi*x)).^2.*(sin(pi*y)).^2.*(sin(pi.*z)).^2);
    end
end