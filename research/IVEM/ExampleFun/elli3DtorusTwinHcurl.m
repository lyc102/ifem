function pde = elli3DtorusTwinHcurl(am,ap,bm,bp,x1,y1,z1,r11,r12,x2,y2,z2,r21,r22)
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
FF1 = @(x,y,z)(((x-x1).^2+(y-y1).^2).^(1/2)-r11).^2+(z-z1).^2-r12^2;
FF2 = @(x,y,z)(((x-x2).^2+(z-z2).^2).^(1/2)-r21).^2+(y-y2).^2-r22^2;
FF = @(x,y,z) ((x-x1).^2+(y-y1).^2).*((x-x2).^2+(z-z2).^2);

scl = 20;

%% interface function
    function u = intf(x,y,z)
        u = min(FF1(x,y,z),FF2(x,y,z));
    end

    function u = intf2(x,y,z)
        u = FF1(x,y,z).*FF2(x,y,z).*FF(x,y,z);
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
        u = DFx(x,y,z)/bm;% + sin(intf(x,y,z))/am;
    end
    function u = um2(x,y,z)
        u = DFy(x,y,z)/bm;% + sin(intf(x,y,z))/am;
    end
    function u = um3(x,y,z)
        u = DFz(x,y,z)/bm + sin(intf2(x,y,z)/scl)/am;
    end
    function u = up1(x,y,z)
        u = DFx(x,y,z)/bp;% + sin(intf(x,y,z))/ap;
    end
    function u = up2(x,y,z)
        u = DFy(x,y,z)/bp;% + sin(intf(x,y,z))/ap;
    end
    function u = up3(x,y,z)
        u = DFz(x,y,z)/bp + sin(intf2(x,y,z)/scl)/ap;
    end
    function u = DFx(x,y,z)
        u = 2*((x - x2).*((x - x1).^2 + (y - y1).^2).*(-r12.^2 + (r11 -...
            ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 + (z - z1).^2).*(-r22.^2 + (y - ...
            y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2) + (x - ...
            x2).*((x - x1).^2 + (y - y1).^2).*(-r12.^2 + (r11 - ...
            ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 + (z - z1).^2).*(-r21 + ...
            ((x - x2).^2 + (z - z2).^2).^(1/2)).*((x - x2).^2 + (z - ...
            z2).^2).^(1/2) + (x - x1).*(-r11 + ...
            ((x - x1).^2 + (y - y1).^2).^(1/2)).*((x - x1).^2 + (y - ...
            y1).^2).^(1/2).*(-r22.^2 + (y - y2).^2 + (r21 - ...
            ((x - x2).^2 + (z - z2).^2).^(1/2)).^2).*((x - x2).^2 + (z - ...
            z2).^2) + (x - x1).*(-r12.^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 +...
            (z - z1).^2).*(-r22.^2 + (y - y2).^2 + (r21 - ...
            ((x - x2).^2 + (z - z2).^2).^(1/2)).^2).*((x - x2).^2 + (z - z2).^2));
    end
    function u = DFy(x,y,z)
        u = 2*(((x - x1).^2 + (y - y1).^2).*(y - ...
            y2).*(-r12.^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 +...
            (z - z1).^2) + (-r11 + ((x - x1).^2 + (y - y1).^2).^(1/2)).*...
            ((x - x1).^2 + (y - y1).^2).^(1/2).*(y - y1).*(-r22.^2 +...
            (y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2) +...
            (y - y1).*(-r12.^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 +...
            (z - z1).^2).*(-r22.^2 + (y - y2).^2 + (r21 - ...
            ((x - x2).^2 + (z - z2).^2).^(1/2)).^2)).*((x - x2).^2 + (z - z2).^2);
    end
    function u = DFz(x,y,z)
        u = 2*((x - x1).^2 + (y - y1).^2).*((z - z1).*(-r22.^2 + (y - y2).^2 +...
            (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2).*((x - x2).^2 +...
            (z - z2).^2) + (-r12.^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 +...
            (z - z1).^2).*(-r22.^2 + (y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2).*...
            (z - z2) + (-r12.^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2+ ...
            (z - z1).^2).*(-r21 + ((x - x2).^2 + (z - z2).^2).^(1/2)).*...
            ((x - x2).^2 + (z - z2).^2).^(1/2).*(z - z2));
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
        u = cos(intf2(x,y,z)/scl).*DFy(x,y,z)/am/scl;
    end
    function u = Dyum(x,y,z)
        u = cos(intf2(x,y,z)/scl).*(-DFx(x,y,z))/am/scl;
    end
    function u = Dzum(x,y,z)
        u = zeros(size(x));
    end
    function u = Dxup(x,y,z)
        u = cos(intf2(x,y,z)/scl).*DFy(x,y,z)/ap/scl;
    end
    function u = Dyup(x,y,z)
        u = cos(intf2(x,y,z)/scl).*(-DFx(x,y,z))/ap/scl;
    end
    function u = Dzup(x,y,z)
        u = zeros(size(x));
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
        u = bm*um1(x,y,z) + cos(intf2(x,y,z)/scl).*DFxz(x,y,z)/scl + ...
            -sin(intf2(x,y,z)/scl).*DFx(x,y,z).*DFz(x,y,z)/scl^2;
    end
    function u = fm2(x,y,z)
        u = bm*um2(x,y,z) +cos(intf2(x,y,z)/scl).*DFyz(x,y,z)/scl + ...
            -sin(intf2(x,y,z)/scl).*DFy(x,y,z).*DFz(x,y,z)/scl^2;
    end
    function u = fm3(x,y,z)
        u = bm*um3(x,y,z) + cos(intf2(x,y,z)/scl).*(-DFxx(x,y,z)-DFyy(x,y,z))/scl + ...
            -sin(intf2(x,y,z)/scl).*(-DFx(x,y,z).*DFx(x,y,z)-DFy(x,y,z).*DFy(x,y,z))/scl^2;
    end
    function u = fp1(x,y,z)
        u = bp*up1(x,y,z) +cos(intf2(x,y,z)/scl).*DFxz(x,y,z)/scl + ...
            -sin(intf2(x,y,z)/scl).*DFx(x,y,z).*DFz(x,y,z)/scl^2;
    end
    function u = fp2(x,y,z)
        u = bp*up2(x,y,z) +cos(intf2(x,y,z)/scl).*DFyz(x,y,z)/scl + ...
            -sin(intf2(x,y,z)/scl).*DFy(x,y,z).*DFz(x,y,z)/scl^2;
    end
    function u = fp3(x,y,z)
        u = bp*up3(x,y,z) +cos(intf2(x,y,z)/scl).*(-DFxx(x,y,z)-DFyy(x,y,z))/scl + ...
            -sin(intf2(x,y,z)/scl).*(-DFx(x,y,z).*DFx(x,y,z)-DFy(x,y,z).*DFy(x,y,z))/scl^2;
    end

    function u = DFxx(x,y,z)
        u = 2.*((x - x2).^2 + (z - z2).^2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2).*...
            ((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) + ((2.*x - 2.*x2).^2.*...
            ((x - x1).^2 + (y - y1).^2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2))/2 +...
            ((2.*x - 2.*x1).^2.*((x - x2).^2 + (z - z2).^2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2))/2 +...
            2.*((x - x1).^2 + (y - y1).^2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2).*...
            ((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) + 2.*(2.*x - 2.*x1).*(2.*x - 2.*x2).*...
            ((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 -...
            r22.^2) - 2.*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*((x - x1).^2 + (y - y1).^2).^(1/2).*...
            ((x - x2).^2 + (z - z2).^2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) -...
            2.*(r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*((x - x1).^2 + (y - y1).^2).*((x - x2).^2 + (z - z2).^2).^(1/2).*...
            ((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2) - (3.*(2.*x - 2.*x1).^2.*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*...
            ((x - x2).^2 + (z - z2).^2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2))./(2.*((x - x1).^2 + (y - y1).^2).^(1/2)) -...
            (3.*(2.*x - 2.*x2).^2.*(r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*((x - x1).^2 + (y - y1).^2).*((z - z1).^2 +...
            (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2))./(2.*((x - x2).^2 + (z - z2).^2).^(1/2)) -...
            2.*(2.*x - 2.*x1).*(2.*x - 2.*x2).*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*((x - x1).^2 + (y - y1).^2).^(1/2).*...
            ((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) - 2.*(2.*x - 2.*x1).*(2.*x - 2.*x2).*(r21 -...
            ((x - x2).^2 + (z - z2).^2).^(1/2)).*((x - x2).^2 + (z - z2).^2).^(1/2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2) +...
            2.*(2.*x - 2.*x1).*(2.*x - 2.*x2).*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*(r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*((x - x1).^2 + (y - y1).^2).^(1/2).*((x - x2).^2 + (z - z2).^2).^(1/2);
    end

    function u = DFyy(x,y,z)
        u = 2.*((x - x2).^2 + (z - z2).^2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2).*((y - y2).^2 +...
            (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) + ((2.*y - 2.*y1).^2.*((x - x2).^2 + (z - z2).^2).*((y - y2).^2 +...
            (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2))/2 + 2.*((x - x1).^2 + (y - y1).^2).*((x - x2).^2 + (z - z2).^2).*...
            ((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2) + 2.*(2.*y - 2.*y1).*(2.*y - 2.*y2).*((x - x2).^2 + (z - z2).^2).*...
            ((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2) - 2.*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*((x - x1).^2 + (y - y1).^2).^(1/2).*...
            ((x - x2).^2 + (z - z2).^2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) - (3.*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*(2.*y - 2.*y1).^2.*...
            ((x - x2).^2 + (z - z2).^2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2))./(2.*((x - x1).^2 + (y - y1).^2).^(1/2)) - 2.*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*...
            (2.*y - 2.*y1).*(2.*y - 2.*y2).*((x - x1).^2 + (y - y1).^2).^(1/2).*((x - x2).^2 + (z - z2).^2);
    end

    function u = DFxz(x,y,z)
        u = ((2.*x - 2.*x2).*(2.*z - 2.*z2).*((x - x1).^2 + (y - y1).^2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2))/2 + (2.*x - 2.*x2).*(2.*z - 2.*z1).*...
            ((x - x1).^2 + (y - y1).^2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) + (2.*x - 2.*x1).*(2.*z - 2.*z1).*((x - x2).^2 + (z - z2).^2).*...
            ((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) + (2.*x - 2.*x1).*(2.*z - 2.*z2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2).*...
            ((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) - (2.*x - 2.*x1).*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*(2.*z - 2.*z2).*((x - x1).^2 + (y - y1).^2).^(1/2).*...
            ((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) - (2.*x - 2.*x1).*(r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*(2.*z - 2.*z2).*((x - x2).^2 + (z - z2).^2).^(1/2).*...
            ((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2) - (2.*x - 2.*x2).*(r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*(2.*z - 2.*z1).*((x - x1).^2 + (y - y1).^2).*...
            ((x - x2).^2 + (z - z2).^2).^(1/2) + (2.*x - 2.*x1).*(r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*(r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*(2.*z - 2.*z2).*((x - x1).^2 + (y - y1).^2).^(1/2).*...
            ((x - x2).^2 + (z - z2).^2).^(1/2) - (3.*(2.*x - 2.*x2).*(r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*(2.*z - 2.*z2).*((x - x1).^2 + (y - y1).^2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2))./(2.*((x - x2).^2 + (z - z2).^2).^(1/2));     
    end

    function u = DFyz(x,y,z)
        u = (2.*y - 2.*y2).*(2.*z - 2.*z2).*((x - x1).^2 + (y - y1).^2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2) + (2.*y - 2.*y1).*(2.*z - 2.*z1).*((x - x2).^2 + (z - z2).^2).*...
            ((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) + (2.*y - 2.*y2).*(2.*z - 2.*z1).*((x - x1).^2 + (y - y1).^2).*((x - x2).^2 + (z - z2).^2) + (2.*y - 2.*y1).*(2.*z - 2.*z2).*...
            ((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) - (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*...
            (2.*y - 2.*y1).*(2.*z - 2.*z2).*((x - x1).^2 + (y - y1).^2).^(1/2).*((y - y2).^2 + (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).^2 - r22.^2) - (2.*y - 2.*y1).*(r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*...
            (2.*z - 2.*z2).*((x - x2).^2 + (z - z2).^2).^(1/2).*((z - z1).^2 + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).^2 - r12.^2) + (r11 - ((x - x1).^2 + (y - y1).^2).^(1/2)).*(2.*y - 2.*y1).*...
            (r21 - ((x - x2).^2 + (z - z2).^2).^(1/2)).*(2.*z - 2.*z2).*((x - x1).^2 + (y - y1).^2).^(1/2).*((x - x2).^2 + (z - z2).^2).^(1/2);
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