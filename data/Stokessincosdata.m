function pde = Stokessincosdata
%
% Created by Ming Wang.
%
nu = 1;
pde = struct('f1',@f1,'f2',@f2,'g',@g,'exactp', @exactp, ...
             'exactu1', @exactu1, 'exactu2', @exactu2,...
             'g_D1', @g_D1, 'g_D2', @g_D2, 'nu', nu);

    function z = f1(p)  
        x = p(:,1); y = p(:,2);
        z = -4*pi^2*(2*cos(2*pi*x)-1).*sin(2*pi*y)+x.^2;

    end

    function z = f2(p)
        x = p(:,1); y = p(:,2);
        z = 4*pi^2*(2*cos(2*pi*y)-1).*sin(2*pi*x);
    end

    function z = g(p)
        z = zeros(size(p,1),1);
    end

    function z = exactu1(p)
        x = p(:,1); y = p(:,2);
        z = (1-cos(2*pi*x)).*sin(2*pi*y); 
    end

    function z = exactu2(p)
        x = p(:,1); y = p(:,2);
        z = -(1-cos(2*pi*y)).*sin(2*pi*x);
    end

    function z = exactp(p)
        x = p(:,1); y = p(:,2);
        z = 1/3*x.^3;
    end

    function z = g_D1(p) % Dirichlet boundary condition
        z = exactu1(p);
    end

    function z = g_D2(p) % Dirichlet boundary condition
        z = exactu2(p);
    end
end