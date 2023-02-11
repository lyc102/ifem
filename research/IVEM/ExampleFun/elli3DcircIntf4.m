function pde = elli3DcircIntf4(am,ap,bm,bp,r1,r2,n1,n2)
%% USAGE: polynomial solution for Poisson equation
%  Last Modified: 02/21/2020 by Xu Zhang

%% PDE Structure
pde = struct('intf',@intf,...
    'exactu1',@exactu1,'exactu2',@exactu2,'exactu3',@exactu3,...
    'um1',@um1,'um2',@um2,'um3',@um3,'up1',@up1,'up2',@up2,'up3',@up3,...
    'Dxu',@Dxu,'Dxum',@Dxum,'Dxup',@Dxup,'Dyu',@Dyu,...
    'Dyum',@Dyum,'Dyup',@Dyup,'Dzu',@Dzu,'Dzum',@Dzum,'Dzup',@Dzup,...
    'f1',@f1,'f2',@f2,'f3',@f3,...
    'fm1',@fm1,'fm2',@fm2,'fm3',@fm3,...
    'fp1',@fp1,'fp2',@fp2,'fp3',@fp3,...
    'A',@A,'Am',@Am,'Ap',@Ap,'one',@one,...
    'B',@B,'Bm',@Bm,'Bp',@Bp);

pde.am = am;
pde.ap = ap;
pde.bm = bm;
pde.bp = bp;
%% interface function
    function u = intf(x,y,z)
        u = (x.^2 + y.^2 + z.^2).^(1/2)/r1-1;
    end

%% exact solution
    function u = exactu1(x,y,z)
        u = um1(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = up1(x(id),y(id),z(id));
    end
    function u = exactu2(x,y,z)
        u = um2(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = up2(x(id),y(id),z(id));
    end
    function u = exactu3(x,y,z)
        u = um3(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = up3(x(id),y(id),z(id));
    end
    function u = um1(x,y,z)
        u = (-y + n1*uker1(x,y,z).*(y-z))/am;
    end
    function u = um2(x,y,z)
        u = (x + n1*uker1(x,y,z).*(z-x))/am;
    end
    function u = um3(x,y,z)
        u = (0 + n1*uker1(x,y,z).*(x-y))/am;
    end
    function u = up1(x,y,z)
        u = (-y + n2*uker1(x,y,z).*uker2(x,y,z).*(y-z))/ap;
    end
    function u = up2(x,y,z)
        u = (x + n2*uker1(x,y,z).*uker2(x,y,z).*(z-x))/ap;
    end
    function u = up3(x,y,z)
        u = (0 + n2*uker1(x,y,z).*uker2(x,y,z).*(x-y))/ap;
    end
% not Hcurl conforming so this one is not correct
    function u = uker1(x,y,z)
        u = r1^2 - (x.^2 + y.^2 + z.^2);
    end
    function u = uker2(x,y,z)
        u = r2^2 - (x.^2 + y.^2 + z.^2);
    end
%% Boundary Function
    function u = gD1(x,y,z)
        u = exactu1(x,y,z);
    end
    function u = gD2(x,y,z)
        u = exactu2(x,y,z);
    end
    function u = gD3(x,y,z)
        u = exactu3(x,y,z);
    end
%% Derivative of the exact solution
    function u = Dxu(x,y,z)
        u = Dxum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dxup(x(id),y(id),z(id));
    end
    function u = Dyu(x,y,z)
        u = Dyum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dyup(x(id),y(id),z(id));
    end
    function u = Dzu(x,y,z)
        u = Dzum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dzup(x(id),y(id),z(id));
    end
    function u = Dxum(x,y,z)
        u = -2*n1/am*(x.*y+x.*z-y.^2-z.^2 + uker1(x,y,z));
    end
    function u = Dyum(x,y,z)
        u = -2*n1/am*(y.*z+x.*y-z.^2-x.^2 + uker1(x,y,z));
    end
    function u = Dzum(x,y,z)
        u = -2*n1/am*(x.*z+y.*z-x.^2-y.^2 + uker1(x,y,z))+2/am;
    end
    function u = Dxup(x,y,z)
        u = -2*n2*uker2(x,y,z)./ap.*(x.*y+x.*z-y.^2-z.^2 + uker1(x,y,z)) - ...
            2*n2*uker1(x,y,z).*(x.*y+x.*z - y.^2-z.^2)/ap;
    end
    function u = Dyup(x,y,z)
        u = -2*n2*uker2(x,y,z)./ap.*(y.*z+x.*y-z.^2-x.^2 + uker1(x,y,z)) - ...
            2*n2*uker1(x,y,z).*(y.*z+x.*y - x.^2-z.^2)/ap;
    end
    function u = Dzup(x,y,z)
        u = -2*n2*uker2(x,y,z)./ap.*(x.*z+y.*z-x.^2-y.^2 + uker1(x,y,z)) - ...
            2*n2*uker1(x,y,z).*(x.*z+y.*z - x.^2-y.^2)/ap+2/ap;
    end

%% right hand side function
    function u = f1(x,y,z)
        u = fm1(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = fp1(x(id),y(id),z(id));
    end
    function u = f2(x,y,z)
        u = fm2(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = fp2(x(id),y(id),z(id));
    end
    function u = f3(x,y,z)
        u = fm3(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = fp3(x(id),y(id),z(id));
    end

    function u = fm1(x,y,z)
        u = -10*n1*(z-y) + bm*um1(x,y,z);
    end
    function u = fm2(x,y,z)
        u = -10*n1*(x-z) + bm*um2(x,y,z);
    end
    function u = fm3(x,y,z)
        u = -10*n1*(y-x) + bm*um3(x,y,z);
    end
    function u = fp1(x,y,z)
        u = -10*n2*uker2(x,y,z).*(z-y) +8*n2*(z.^3-y.^3+z.*y.^2+z.*x.^2-y.*x.^2-y.*z.^2) -...
             10*n2*uker1(x,y,z).*(z-y) + bp*up1(x,y,z);
    end
    function u = fp2(x,y,z)
        u = -10*n2*uker2(x,y,z).*(x-z) +8*n2*(x.^3-z.^3+x.*z.^2+x.*y.^2-z.*y.^2-z.*x.^2) -...
             10*n2*uker1(x,y,z).*(x-z) + bp*up2(x,y,z);
    end
    function u = fp3(x,y,z)
        u = -10*n2*uker2(x,y,z).*(y-x) +8*n2*(y.^3-x.^3+y.*x.^2+y.*z.^2-x.*z.^2-x.*y.^2) -...
             10*n2*uker1(x,y,z).*(y-x) + bp*up3(x,y,z);
    end

%% Diffusion coefficient function
    function u = A(x,y,z)
        u = Am(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Ap(x(id),y(id),z(id));
    end
    function u = Am(x,y,z)
        u = am*ones(size(x));
    end
    function u = Ap(x,y,z)
        u = ap*ones(size(x));
    end

%% Mass coefficient function
    function u = B(x,y,z)
        u = Bm(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Bp(x(id),y(id),z(id));
    end
    function u = Bm(x,y,z)
        u = bm*ones(size(x));
    end
    function u = Bp(x,y,z)
        u = bp*ones(size(x));
    end

%% Other function
    function u = one(x,y,z)
        u = ones(size(x));
    end
end