function pde = elli3DlinIntf(bm,bp,cx,cy,cz,rx,ry,rz,a,cm,cp)
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
        u = ((x-cx).*rx + (y-cy).*ry + (z-cz).*rz)/(rx^2+ry^2+rz^2)^(1/2);
    end

%% exact solution
    function u = exactu(x,y,z)
        u = um(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = up(x(id),y(id),z(id));
    end
    function u = um(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = cm*r.^a/bm;
    end
    function u = up(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = cp*r.^a/bp;
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
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = cm*a*r.^(a-1)/bm*rx;
    end
    function u = Dxup(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = cp*a*r.^(a-1)/bp*rx;
    end

    function u = Dyu(x,y,z)
        u = Dyum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dyup(x(id),y(id),z(id));
    end
    function u = Dyum(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = cm*a*r.^(a-1)/bm*ry;
    end
    function u = Dyup(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = cp*a*r.^(a-1)/bp*ry;
    end

    function u = Dzu(x,y,z)
        u = Dzum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dzup(x(id),y(id),z(id));
    end
    function u = Dzum(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = cm*a*r.^(a-1)/bm*rz;
    end
    function u = Dzup(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = cp*a*r.^(a-1)/bp*rz;
    end

%% right hand side function
    function u = f(x,y,z)
        u = fm(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = fp(x(id),y(id),z(id));
    end
    function u = fm(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = -a*(a-1)*r.^(a-2).*(rx^2 + ry^2 + rz^2);
    end
    function u = fp(x,y,z)
        r = (x-cx).*rx + (y-cy).*ry +(z-cz).*rz;
        u = -a*(a-1)*r.^(a-2).*(rx^2 + ry^2 + rz^2);
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