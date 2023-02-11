function pde = elli3DtorusTwin(bm,bp,x1,y1,z1,r11,r12,x2,y2,z2,r21,r22)
%% USAGE: polynomial solution for Poisson equation
%  Last Modified: 02/21/2020 by Xu Zhang

%% PDE Structure
pde = struct('intf',@intf,'f',@f,'fm',@fm,'fp',@fp,'exactu',@exactu,...
    'um',@um,'up',@up,'Dxu',@Dxu,'Dxum',@Dxum,'Dxup',@Dxup,'Dyu',@Dyu,...
    'Dyum',@Dyum,'Dyup',@Dyup,'Dzu',@Dzu,'Dzum',@Dzum,'Dzup',@Dzup,...
    'A',@A,'Am',@Am,'Ap',@Ap,'one',@one,'gD',@gD);

pde.bm = bm;
pde.bp = bp;
bbm =1; bbp = 1;
FF1 = @(x,y,z)(((x-x1).^2+(y-y1).^2).^(1/2)-r11).^2+(z-z1).^2-r12^2;
FF2 = @(x,y,z)(((x-x2).^2+(z-z2).^2).^(1/2)-r21).^2+(y-y2).^2-r22^2;

%% interface function
    function u = intf(x,y,z)
        u = min(FF1(x,y,z),FF2(x,y,z));
    end

    function u = intf2(x,y,z)
        u = FF1(x,y,z).*FF2(x,y,z);
    end

%% exact solution
    function u = exactu(x,y,z)
        u = um(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = up(x(id),y(id),z(id));
    end
    function u = um(x,y,z)
        u = ones(size(x));
    end
    function u = up(x,y,z)
        u = cos(intf2(x,y,z)/bbp);
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
        u = zeros(size(x));
    end
    function u = Dxup(x,y,z)
        u = -sin(intf2(x,y,z)/bbp).*DFx(x,y,z)/bbp;
    end

    function u = Dyu(x,y,z)
        u = Dyum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dyup(x(id),y(id),z(id));
    end
    function u = Dyum(x,y,z)
        u =  zeros(size(x));
    end
    function u = Dyup(x,y,z)
        u = -sin(intf2(x,y,z)/bbp).*DFy(x,y,z)/bbp;
    end

    function u = Dzu(x,y,z)
        u = Dzum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Dzup(x(id),y(id),z(id));
    end
    function u = Dzum(x,y,z)
        u = zeros(size(x));
    end
    function u = Dzup(x,y,z)
        u = -sin(intf2(x,y,z)/bbp).*DFz(x,y,z)/bbp;
    end

    function u = DFx(x,y,z)
       u = 2*( ((x-x1).^2+(y-y1).^2).^(1/2)-r11 ).*((x-x1).^2+(y-y1).^2).^(-1/2).*...
           (x-x1).*FF2(x,y,z)+...
           2*( ((x-x2).^2+(z-z2).^2).^(1/2)-r21 ).*((x-x2).^2+(z-z2).^2).^(-1/2).*...
           (x-x2).*FF1(x,y,z);
    end
    function u = DFy(x,y,z)
        u = 2*( ((x-x1).^2+(y-y1).^2).^(1/2)-r11 ).*((x-x1).^2+(y-y1).^2).^(-1/2).*...
            (y-y1).*FF2(x,y,z)+...
            2*(y-y2).*FF1(x,y,z);
    end
    function u = DFz(x,y,z)
        u = 2*( ((x-x2).^2+(z-z2).^2).^(1/2)-r21 ).*((x-x2).^2+(z-z2).^2).^(-1/2).*...
            (z-z2).*FF1(x,y,z)+...
            2*(z-z1).*FF2(x,y,z);
    end

%% right hand side function
    function u = f(x,y,z)
        u = fm(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = fp(x(id),y(id),z(id));
    end
    function u = fm(x,y,z)
        u = zeros(size(x));
    end
    function u = fp(x,y,z)
        u = (sin(intf2(x,y,z)/bbp).*(DFxx(x,y,z)+DFyy(x,y,z)+DFzz(x,y,z))+...
            cos(intf2(x,y,z)/bbp).*(DFx(x,y,z).^2 + DFy(x,y,z).^2 + DFz(x,y,z).^2)/bbp)*bp;
    end

    function u = DFxx(x,y,z)
        u = 2*( (x-x1).*((x-x1).^2+(y-y1).^2).^(-1).*(x-x1)+...
            -(((x-x1).^2+(y-y1).^2).^(1/2)-r11).*((x-x1).^2+(y-y1).^2).^(-3/2).*(x-x1).^2+...
            (((x-x1).^2+(y-y1).^2).^(1/2)-r11).*((x-x1).^2+(y-y1).^2).^(-1/2) ).*FF2(x,y,z)+...
            2*( (x-x2).*((x-x2).^2+(z-z2).^2).^(-1).*(x-x2)+...
            -(((x-x2).^2+(z-z2).^2).^(1/2)-r21).*((x-x2).^2+(z-z2).^2).^(-3/2).*(x-x2).^2+...
            (((x-x2).^2+(z-z2).^2).^(1/2)-r21).*((x-x2).^2+(z-z2).^2).^(-1/2) ).*FF1(x,y,z) +...
            8*( ((x-x1).^2+(y-y1).^2).^(1/2)-r11 ).*((x-x1).^2+(y-y1).^2).^(-1/2).*(x-x1).*...
            ( ((x-x2).^2+(z-z2).^2).^(1/2)-r21 ).*((x-x2).^2+(z-z2).^2).^(-1/2).*(x-x2);
    end
    function u = DFyy(x,y,z)
        u = 2*( (y-y1).*((x-x1).^2+(y-y1).^2).^(-1).*(y-y1)+...
            -(((x-x1).^2+(y-y1).^2).^(1/2)-r11).*((x-x1).^2+(y-y1).^2).^(-3/2).*(y-y1).^2+...
            (((x-x1).^2+(y-y1).^2).^(1/2)-r11).*((x-x1).^2+(y-y1).^2).^(-1/2) ).*FF2(x,y,z)+...
            2*FF1(x,y,z) +...
            8*( ((x-x1).^2+(y-y1).^2).^(1/2)-r11 ).*((x-x1).^2+(y-y1).^2).^(-1/2).*(y-y1).*(y-y2);
    end
    function u = DFzz(x,y,z)
        u = 2*( (z-z2).*((x-x2).^2+(z-z2).^2).^(-1).*(z-z2)+...
            -(((x-x2).^2+(z-z2).^2).^(1/2)-r21).*((x-x2).^2+(z-z2).^2).^(-3/2).*(z-z2).^2+...
            (((x-x2).^2+(z-z2).^2).^(1/2)-r21).*((x-x2).^2+(z-z2).^2).^(-1/2) ).*FF1(x,y,z)+...
            2*FF2(x,y,z) +...
            8*( ((x-x2).^2+(z-z2).^2).^(1/2)-r21 ).*((x-x2).^2+(z-z2).^2).^(-1/2).*(z-z2).*(z-z1);
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