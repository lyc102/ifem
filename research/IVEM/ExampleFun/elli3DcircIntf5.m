function pde = elli3DcircIntf5(bm,bp,r,x0,y0,z0)
%% USAGE: polynomial solution for Poisson equation
%  Last Modified: 02/21/2020 by Xu Zhang

%% PDE Structure
pde = struct('intf',@intf,'f',@f,'fm',@fm,'fp',@fp,'exactu',@exactu,...
    'um',@um,'up',@up,'Dxu',@Dxu,'Dxum',@Dxum,'Dxup',@Dxup,'Dyu',@Dyu,...
    'Dyum',@Dyum,'Dyup',@Dyup,'Dzu',@Dzu,'Dzum',@Dzum,'Dzup',@Dzup,...
    'A',@A,'Am',@Am,'Ap',@Ap,'one',@one,'gD',@gD);

pde.bm = bm;
pde.bp = bp;
%% interface function
    function u = intf(x,y,z)
        u = ((x-x0).^2/r^2 + (y-y0).^2/r^2 + (z-z0).^2/r^2).^(1/2)-1;
    end

%% exact solution
    function u = exactu(x,y,z)
        u = um(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = up(x(id),y(id),z(id));
    end
    function u = um(x,y,z)
        u = exp((x.^2+y.^2+z.^2-r^2)/bm);
    end
    function u = up(x,y,z)
        u = sin((x.^2+y.^2+z.^2-r^2)/bp)+1;
    end
%% Boundary Function
    function u = gD(x,y,z)
        u = exactu(x,y,z);
    end
%% Derivative of the exact solution
    function u = Dxu(x,y,z)
        u = Dxum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dxup(x(id),y(id),z(id));
    end
    function u = Dxum(x,y,z)
        u = exp((x.^2+y.^2+z.^2-r^2)/bm).*(2*x)/bm;
    end
    function u = Dxup(x,y,z)
        u = cos((x.^2+y.^2+z.^2-r^2)/bp).*(2*x)/bp;
    end

    function u = Dyu(x,y,z)
        u = Dyum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dyup(x(id),y(id),z(id));
    end
    function u = Dyum(x,y,z)
        u = exp((x.^2+y.^2+z.^2-r^2)/bm).*(2*y)/bm;
    end
    function u = Dyup(x,y,z)
        u = cos((x.^2+y.^2+z.^2-r^2)/bp).*(2*y)/bp;
    end

    function u = Dzu(x,y,z)
        u = Dzum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dzup(x(id),y(id),z(id));
    end
    function u = Dzum(x,y,z)
        u = exp((x.^2+y.^2+z.^2-r^2)/bm).*(2*z)/bm;
    end
    function u = Dzup(x,y,z)
        u = cos((x.^2+y.^2+z.^2-r^2)/bp).*(2*z)/bp;
    end

%% right hand side function
    function u = f(x,y,z)
        u = fm(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = fp(x(id),y(id),z(id));
    end
    function u = fm(x,y,z)
        u = -6*exp((x.^2+y.^2+z.^2-r^2)/bm) +...
            -4*exp((x.^2+y.^2+z.^2-r^2)/bm).*(x.^2+y.^2+z.^2)/bm;
    end
    function u = fp(x,y,z)
        u = -6*cos((x.^2+y.^2+z.^2-r^2)/bp) + ...
            4*sin((x.^2+y.^2+z.^2-r^2)/bp).*(x.^2+y.^2+z.^2)/bp;
    end

%% Diffusion coefficient function
    function u = A(x,y,z)
        u = Am(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Ap(x(id),y(id),z(id));
    end
    function u = Am(x,y,z)
        u = bm*ones(size(x));
    end
    function u = Ap(x,y,z)
        u = bp*ones(size(x));
    end

%% Other function
    function u = one(x,y,z)
        u = ones(size(x));
    end
end