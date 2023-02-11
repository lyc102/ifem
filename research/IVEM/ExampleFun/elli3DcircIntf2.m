function pde = elli3DcircIntf2(am,ap,bm,bp,r,x0,y0,z0,a11,a12,a)
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
        u = ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(1/2)/r-1;
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
        u1 = uker1(x,y,z).*(x-x0);
        u2 = uker2(x,y,z);
        u = (u1 + u2)/bm;
    end
    function u = um2(x,y,z)
        u1 = uker1(x,y,z).*(y-y0);
        u2 = uker2(x,y,z);
        u = (u1 + u2)/bm;
    end
    function u = um3(x,y,z)
        u1 = uker1(x,y,z).*(z-z0);
        u2 = uker2(x,y,z);
        u = (u1 + u2)/bm;
    end
    function u = up1(x,y,z)
        u1 = uker1(x,y,z).*(x-x0);
        u2 = uker2(x,y,z);
        u = (u1 + u2)/bp;
    end
    function u = up2(x,y,z)
        u1 = uker1(x,y,z).*(y-y0);
        u2 = uker2(x,y,z);
        u = (u1 + u2)/bp;
    end
    function u = up3(x,y,z)
        u1 = uker1(x,y,z).*(z-z0);
        u2 = uker2(x,y,z);
        u = (u1 + u2)/bp;
    end
    function u = uker1(x,y,z)
        u = nthroot((x-x0).^2 + (y-y0).^2 + (z-z0).^2 - r^2,a11);
        u = u.^a12;
    end
    function u = uker2(x,y,z)
        u = ((x-x0).^2 + (y-y0).^2 + (z-z0).^2 - r^2).^a;
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
        u = (Duker(x,y,z).*(y-y0)*2 - Duker(x,y,z).*(z-z0)*2)/bm;
    end
    function u = Dyum(x,y,z)
        u = (Duker(x,y,z).*(z-z0)*2 - Duker(x,y,z).*(x-x0)*2)/bm;
    end
    function u = Dzum(x,y,z)
        u = (Duker(x,y,z).*(x-x0)*2 - Duker(x,y,z).*(y-y0)*2)/bm;
    end
    function u = Dxup(x,y,z)
        u = (Duker(x,y,z).*(y-y0)*2 - Duker(x,y,z).*(z-z0)*2)/bp;
    end
    function u = Dyup(x,y,z)
        u = (Duker(x,y,z).*(z-z0)*2 - Duker(x,y,z).*(x-x0)*2)/bp;
    end
    function u = Dzup(x,y,z)
        u = (Duker(x,y,z).*(x-x0)*2 - Duker(x,y,z).*(y-y0)*2)/bp;
    end

    function u = Duker(x,y,z)
        u = a*((x-x0).^2 + (y-y0).^2 + (z-z0).^2 - r^2).^(a-1);
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
        xh = x-x0; yh = y-y0; zh = z-z0;
        u = (-4*Duker(x,y,z) + DDuker(x,y,z).*(xh.*(yh+zh)-yh.^2-zh.^2)) + bm*um1(x,y,z);
    end
    function u = fm2(x,y,z)
        xh = x-x0; yh = y-y0; zh = z-z0;
        u = (-4*Duker(x,y,z) + DDuker(x,y,z).*(yh.*(xh+zh)-xh.^2-zh.^2)) + bm*um2(x,y,z);
    end
    function u = fm3(x,y,z)
        xh = x-x0; yh = y-y0; zh = z-z0;
        u = (-4*Duker(x,y,z) + DDuker(x,y,z).*(zh.*(xh+yh)-xh.^2-yh.^2)) + bm*um3(x,y,z);
    end
    function u = fp1(x,y,z)
        xh = x-x0; yh = y-y0; zh = z-z0;
        u = (-4*Duker(x,y,z) + DDuker(x,y,z).*(xh.*(yh+zh)-yh.^2-zh.^2)) + bp*up1(x,y,z);
    end
    function u = fp2(x,y,z)
        xh = x-x0; yh = y-y0; zh = z-z0;
        u = (-4*Duker(x,y,z) + DDuker(x,y,z).*(yh.*(xh+zh)-xh.^2-zh.^2)) + bp*up2(x,y,z);
    end
    function u = fp3(x,y,z)
        xh = x-x0; yh = y-y0; zh = z-z0;
        u = (-4*Duker(x,y,z) + DDuker(x,y,z).*(zh.*(xh+yh)-xh.^2-yh.^2)) + bp*up3(x,y,z);
    end

    function u = DDuker(x,y,z)
        u = a*(a-1)*((x-x0).^2 + (y-y0).^2 + (z-z0).^2 - r^2).^(a-2);
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