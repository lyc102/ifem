function pde = elli3DlinIntf2(bm,bp,delta)
%% USAGE: polynomial solution for Poisson equation
%  Last Modified: 12/15/2020 by GRC

%% PDE Structure
pde = struct('intf',@intf,'f',@f,'fm',@fm,'fp',@fp,'exactu',@exactu,...
    'um',@um,'up',@up,'Dxu',@Dxu,'Dxum',@Dxum,'Dxup',@Dxup,'Dyu',@Dyu,...
    'Dyum',@Dyum,'Dyup',@Dyup,'Dzu',@Dzu,'Dzum',@Dzum,'Dzup',@Dzup,...
    'A',@A,'Am',@Am,'Ap',@Ap,'one',@one,'gD',@gD);

pde.bm = bm;
pde.bp = bp;
%% interface function
    function u = intf(x,y,z)
        u = (x-delta);
    end

%% exact solution
    function u = exactu(x,y,z)
        u = um(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = up(x(id),y(id),z(id));
    end
    function u = um(x,y,z)
        r = (1+delta)*sin((x-delta)*pi/(1+delta)).*sin(pi*y).*sin(pi*z);
        u = r/bm;
    end
    function u = up(x,y,z)
        r = (1-delta)*sin((x-delta)*pi/(1-delta)).*sin(pi*y).*sin(pi*z);
        u = r/bp;
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
        r = cos((x-delta)*pi/(1+delta)).*sin(pi*y).*sin(pi*z)*pi;
        u = r/bm;
    end
    function u = Dxup(x,y,z)
        r = cos((x-delta)*pi/(1-delta)).*sin(pi*y).*sin(pi*z)*pi;
        u = r/bp;
    end

    function u = Dyu(x,y,z)
        u = Dyum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dyup(x(id),y(id),z(id));
    end
    function u = Dyum(x,y,z)
        r = (1+delta)*sin((x-delta)*pi/(1+delta)).*cos(pi*y).*sin(pi*z)*pi;
        u = r/bm;
    end
    function u = Dyup(x,y,z)
        r = (1-delta)*sin((x-delta)*pi/(1-delta)).*cos(pi*y).*sin(pi*z)*pi;
        u = r/bp;
    end

    function u = Dzu(x,y,z)
        u = Dzum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dzup(x(id),y(id),z(id));
    end
    function u = Dzum(x,y,z)
        r = (1+delta)*sin((x-delta)*pi/(1+delta)).*sin(pi*y).*cos(pi*z)*pi;
        u = r/bm;
    end
    function u = Dzup(x,y,z)
        r = (1-delta)*sin((x-delta)*pi/(1-delta)).*sin(pi*y).*cos(pi*z)*pi;
        u = r/bp;
    end

%% right hand side function
    function u = f(x,y,z)
        u = fm(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = fp(x(id),y(id),z(id));
    end
    function u = fm(x,y,z)
        r = sin((x-delta)*pi/(1+delta)).*sin(pi*y).*sin(pi*z)*pi^2;
        u = r*(1/(1+delta)+2*(1+delta));
    end
    function u = fp(x,y,z)
        r = sin((x-delta)*pi/(1-delta)).*sin(pi*y).*sin(pi*z)*pi^2;
        u = r*(1/(1-delta)+2*(1-delta));
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