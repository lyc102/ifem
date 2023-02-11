function pde = elli3DorthocircIntf(bm,bp,rx,ry,rz)
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
        u = log(F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1);
    end

%% exact solution
    function u = exactu(x,y,z)
        u = um(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = up(x(id),y(id),z(id));
    end
    function u = um(x,y,z)
        u = intf(x,y,z)/bm;
    end
    function u = up(x,y,z)
        u = intf(x,y,z)/bp;
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
        fx = (F1x(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2x(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3x(x,y,z)) - Dxr(x,y,z);
        u = (F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1).^(-1).*fx/bm;
    end
    function u = Dxup(x,y,z)
        fx = (F1x(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2x(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3x(x,y,z)) - Dxr(x,y,z);
        u = (F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1).^(-1).*fx/bp;
    end

    function u = Dyu(x,y,z)
        u = Dyum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dyup(x(id),y(id),z(id));
    end
    function u = Dyum(x,y,z)
        fy = (F1y(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2y(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3y(x,y,z)) - Dyr(x,y,z);
        u = (F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1).^(-1).*fy/bm;
    end
    function u = Dyup(x,y,z)
        fy = (F1y(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2y(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3y(x,y,z)) - Dyr(x,y,z);
        u = (F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1).^(-1).*fy/bp;
    end

    function u = Dzu(x,y,z)
        u = Dzum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dzup(x(id),y(id),z(id));
    end
    function u = Dzum(x,y,z)
        fz = (F1z(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2z(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3z(x,y,z)) - Dzr(x,y,z);
        u = (F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1).^(-1).*fz/bm;
    end
    function u = Dzup(x,y,z)
        fz = (F1z(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2z(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3z(x,y,z)) - Dzr(x,y,z);
        u = (F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1).^(-1).*fz/bp;
    end

%% right hand side function
    function u = f(x,y,z)
        u = fm(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = fp(x(id),y(id),z(id));
    end
    function u = fm(x,y,z)
        r_xx = ry^2*(rz*2);
        r_yy = ry^2*(rz*2);
        r_zz = ry^2*(rz*2);
        
        fx = (F1x(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2x(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3x(x,y,z)) - Dxr(x,y,z);
        fy = (F1y(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2y(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3y(x,y,z)) - Dyr(x,y,z);
        fz = (F1z(x,y,z).*F2(x,y,z).*F3(x,y,z)+F1(x,y,z).*F2z(x,y,z).*F3(x,y,z)+...
            F1(x,y,z).*F2(x,y,z).*F3z(x,y,z)) - Dzr(x,y,z);
        
        fxx = F1xx(x,y,z).*F2(x,y,z).*F3(x,y,z) + F1x(x,y,z).*F2x(x,y,z).*F3(x,y,z) + F1x(x,y,z).*F2(x,y,z).*F3x(x,y,z) +...
            F1(x,y,z).*F2xx(x,y,z).*F3(x,y,z) + F1x(x,y,z).*F2x(x,y,z).*F3(x,y,z) + F1(x,y,z).*F2x(x,y,z).*F3x(x,y,z) +...
            F1(x,y,z).*F2(x,y,z).*F3xx(x,y,z) + F1(x,y,z).*F2x(x,y,z).*F3x(x,y,z) + F1x(x,y,z).*F2(x,y,z).*F3x(x,y,z) - r_xx;
        
        fyy = F1yy(x,y,z).*F2(x,y,z).*F3(x,y,z) + F1y(x,y,z).*F2y(x,y,z).*F3(x,y,z) + F1y(x,y,z).*F2(x,y,z).*F3y(x,y,z) +...
            F1(x,y,z).*F2yy(x,y,z).*F3(x,y,z) + F1y(x,y,z).*F2y(x,y,z).*F3(x,y,z) + F1(x,y,z).*F2y(x,y,z).*F3y(x,y,z) +...
            F1(x,y,z).*F2(x,y,z).*F3yy(x,y,z) + F1(x,y,z).*F2y(x,y,z).*F3y(x,y,z) + F1y(x,y,z).*F2(x,y,z).*F3y(x,y,z) - r_yy;
        
        fzz = F1zz(x,y,z).*F2(x,y,z).*F3(x,y,z) + F1z(x,y,z).*F2z(x,y,z).*F3(x,y,z) + F1z(x,y,z).*F2(x,y,z).*F3z(x,y,z) +...
            F1(x,y,z).*F2zz(x,y,z).*F3(x,y,z) + F1z(x,y,z).*F2z(x,y,z).*F3(x,y,z) + F1(x,y,z).*F2z(x,y,z).*F3z(x,y,z) +...
            F1(x,y,z).*F2(x,y,z).*F3zz(x,y,z) + F1(x,y,z).*F2z(x,y,z).*F3z(x,y,z) + F1z(x,y,z).*F2(x,y,z).*F3z(x,y,z) - r_zz;
        
        u = -(F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1).^(-1).*(fxx+fyy+fzz) +...
            (F1(x,y,z).*F2(x,y,z).*F3(x,y,z)-r(x,y,z)+1).^(-2).*(fx.^2+fy.^2+fz.^2);
    end
    function u = fp(x,y,z)
        u = fm(x,y,z);
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
    function u = r(x,y,z)
        u = ry^2*(1+rz*(x.^2 + y.^2 + z.^2));
    end
    function u = Dxr(x,y,z)
        u = ry^2*(rz*2*x);
    end
    function u = Dyr(x,y,z)
        u = ry^2*(rz*2*y);
    end
    function u = Dzr(x,y,z)
        u = ry^2*(rz*2*z);
    end
    function u = F1(x,y,z)
        u = (x.^2 + y.^2 -rx^2).^2 + z.^2;
    end
    function u = F2(x,y,z)
        u = (x.^2 + z.^2 -rx^2).^2 + y.^2;
    end
    function u = F3(x,y,z)
        u = (y.^2 + z.^2 -rx^2).^2 + x.^2;
    end
    function u = F1x(x,y,z)
        u = 2*(x.^2 + y.^2 -rx^2).*(2*x);
    end
    function u = F1xx(x,y,z)
        u = 2*(x.^2 + y.^2 -rx^2).*(2) + 8*x.^2;
    end
    function u = F2x(x,y,z)
        u = 2*(x.^2 + z.^2 -rx^2).*(2*x);
    end
    function u = F2xx(x,y,z)
        u = 2*(x.^2 + z.^2 -rx^2).*(2) + 8*x.^2;
    end
    function u = F3x(x,y,z)
        u = 2*x;
    end
    function u = F3xx(x,y,z)
        u = 2*ones(size(x));
    end
    function u = F1y(x,y,z)
        u = 2*(x.^2 + y.^2 -rx^2).*(2*y);
    end
    function u = F1yy(x,y,z)
        u = 2*(x.^2 + y.^2 -rx^2).*(2) + 8*y.^2;
    end
    function u = F2y(x,y,z)
        u = 2*y;
    end
    function u = F2yy(x,y,z)
        u = 2*ones(size(y));
    end
    function u = F3y(x,y,z)
        u = 2*(y.^2 + z.^2 -rx^2).*(2*y);
    end
    function u = F3yy(x,y,z)
        u = 2*(y.^2 + z.^2 -rx^2)*(2) + 8*y.^2;
    end
    function u = F1z(x,y,z)
        u = 2*z;
    end
    function u = F1zz(x,y,z)
        u = 2*ones(size(z));
    end
    function u = F2z(x,y,z)
        u = 2*(x.^2 + z.^2 -rx^2).*(2*z);
    end
    function u = F2zz(x,y,z)
        u = 2*(x.^2 + z.^2 -rx^2).*(2) + 8*z.^2;
    end
    function u = F3z(x,y,z)
        u = 2*(y.^2 + z.^2 -rx^2).*(2*z);
    end
    function u = F3zz(x,y,z)
        u = 2*(y.^2 + z.^2 -rx^2).*(2) + 8*z.^2;
    end
end