function pde = Stokesdata0

%
nu = 1;
pde = struct('f',@f,'g',@g,'exactp', @exactp, ...
             'exactu', @exactu,'g_D', @g_D, 'nu', nu);

    function z = f(p)  
        x = p(:,1); y = p(:,2);
        z(:,1) = -4*pi^2*(2*cos(2*pi*x)-1).*sin(2*pi*y)+x.^2;
        z(:,2) = 4*pi^2*(2*cos(2*pi*y)-1).*sin(2*pi*x);
    end

    function z = g(p)
        z = zeros(size(p,1),1);
    end

    function z = exactu(p)
        x = p(:,1); y = p(:,2);
        z(:,1) = (1-cos(2*pi*x)).*sin(2*pi*y); 
        z(:,2) = -(1-cos(2*pi*y)).*sin(2*pi*x);
    end

    function z = exactp(p)
        x = p(:,1); % y = p(:,2);
        z = 1/3*x.^3;
    end

    function z = g_D(p) % Dirichlet boundary condition
        z = exactu(p);
    end
end