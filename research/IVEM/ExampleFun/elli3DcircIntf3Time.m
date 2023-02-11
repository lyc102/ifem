function pde = elli3DcircIntf3Time(mum,mup,sigm,sigp,epsm,epsp,omega,r1,r2,n1,n2)
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
    'Mu',@Mu,'Mum',@Mum,'Mup',@Mup,'one',@one,...
    'Epslon',@Epslon,'Epslonm',@Epslonm,'Epslonp',@Epslonp,...
    'Sig',@Sig,'Sigm',@Sigm,'Sigp',@Sigp);

pde.mum = mum;
pde.mup = mup;
pde.sigm = sigm;
pde.sigp = sigp;
pde.epsm = epsm;
pde.epsp = epsp;
%% interface function
    function u = intf(x,y,z)
        u = (x.^2 + y.^2 + z.^2).^(1/2)/r1-1;
    end

%% exact solution
    function u = exactu1(x,y,z,t)
        u = um1(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = up1(x(id),y(id),z(id),t);
        %u = t*ones(size(x));
    end
    function u = exactu2(x,y,z,t)
        u = um2(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = up2(x(id),y(id),z(id),t);
        %u = t*ones(size(x));
    end
    function u = exactu3(x,y,z,t)
        u = um3(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = up3(x(id),y(id),z(id),t);
        %u = t*ones(size(x));
    end
    function u = um1(x,y,z,t)
        u = (x + n1*uker1(x,y,z).*(y-z))/mum*sin(2*pi*omega*t);
        %u = t*ones(size(x));
    end
    function u = um2(x,y,z,t)
        u = (y + n1*uker1(x,y,z).*(z-x))/mum*sin(2*pi*omega*t);
        %u = t*ones(size(x));
    end
    function u = um3(x,y,z,t)
        u = (z + n1*uker1(x,y,z).*(x-y))/mum*sin(2*pi*omega*t);
        %u = t*ones(size(x));
    end
    function u = up1(x,y,z,t)
        u = (x + n2*uker1(x,y,z).*uker2(x,y,z).*(y-z))/mup*sin(2*pi*omega*t);
        %u = t*ones(size(x));
    end
    function u = up2(x,y,z,t)
        u = (y + n2*uker1(x,y,z).*uker2(x,y,z).*(z-x))/mup*sin(2*pi*omega*t);
        %u = t*ones(size(x));
    end
    function u = up3(x,y,z,t)
        u = (z + n2*uker1(x,y,z).*uker2(x,y,z).*(x-y))/mup*sin(2*pi*omega*t);
        %u = t*ones(size(x));
    end
    function u = uker1(x,y,z)
        u = r1^2 - (x.^2 + y.^2 + z.^2);
    end
    function u = uker2(x,y,z)
        u = r2^2 - (x.^2 + y.^2 + z.^2);
    end
%% Boundary Function
    function u = gD1(x,y,z,t)
        u = exactu1(x,y,z,t);
    end
    function u = gD2(x,y,z,t)
        u = exactu2(x,y,z,t);
    end
    function u = gD3(x,y,z,t)
        u = exactu3(x,y,z,t);
    end
%% Derivative of the exact solution
    function u = Dxu(x,y,z,t)
        u = Dxum(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = Dxup(x(id),y(id),z(id),t);
    end
    function u = Dyu(x,y,z,t)
        u = Dyum(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = Dyup(x(id),y(id),z(id),t);
    end
    function u = Dzu(x,y,z,t)
        u = Dzum(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = Dzup(x(id),y(id),z(id),t);
    end
    function u = Dxum(x,y,z,t)
        u = -2*n1/mum*(x.*y+x.*z-y.^2-z.^2 + uker1(x,y,z))*sin(2*pi*omega*t);
    end
    function u = Dyum(x,y,z,t)
        u = -2*n1/mum*(y.*z+x.*y-z.^2-x.^2 + uker1(x,y,z))*sin(2*pi*omega*t);
    end
    function u = Dzum(x,y,z,t)
        u = -2*n1/mum*(x.*z+y.*z-x.^2-y.^2 + uker1(x,y,z))*sin(2*pi*omega*t);
    end
    function u = Dxup(x,y,z,t)
        u = (-2*n2*uker2(x,y,z)./mup.*(x.*y+x.*z-y.^2-z.^2 + uker1(x,y,z)) - ...
            2*n2*uker1(x,y,z).*(x.*y+x.*z - y.^2-z.^2)/mup)*sin(2*pi*omega*t);
    end
    function u = Dyup(x,y,z,t)
        u = (-2*n2*uker2(x,y,z)./mup.*(y.*z+x.*y-z.^2-x.^2 + uker1(x,y,z)) - ...
            2*n2*uker1(x,y,z).*(y.*z+x.*y - x.^2-z.^2)/mup)*sin(2*pi*omega*t);
    end
    function u = Dzup(x,y,z,t)
        u = (-2*n2*uker2(x,y,z)./mup.*(x.*z+y.*z-x.^2-y.^2 + uker1(x,y,z)) - ...
            2*n2*uker1(x,y,z).*(x.*z+y.*z - x.^2-y.^2)/mup)*sin(2*pi*omega*t);
    end

%% right hand side function
    function u = f1(x,y,z,t)
        u = fm1(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = fp1(x(id),y(id),z(id),t);
    end
    function u = f2(x,y,z,t)
        u = fm2(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = fp2(x(id),y(id),z(id),t);
    end
    function u = f3(x,y,z,t)
        u = fm3(x,y,z,t);
        id = intf(x,y,z) > 0;
        u(id) = fp3(x(id),y(id),z(id),t);
    end

    function u = fm1(x,y,z,t)
        u = -10*n1*(z-y)*sin(2*pi*omega*t) - 4*pi^2*omega^2*epsm*um1(x,y,z,t) +...
            2*pi*omega*sigm*(x + n1*uker1(x,y,z).*(y-z))/mum*cos(2*pi*omega*t);
        %u = ones(size(x));
    end
    function u = fm2(x,y,z,t)
        u = -10*n1*(x-z)*sin(2*pi*omega*t) - 4*pi^2*omega^2*epsm*um2(x,y,z,t) +...
            2*pi*omega*sigm*(y + n1*uker1(x,y,z).*(z-x))/mum*cos(2*pi*omega*t);
        %u = ones(size(x));
    end
    function u = fm3(x,y,z,t)
        u = -10*n1*(y-x)*sin(2*pi*omega*t) - 4*pi^2*omega^2*epsm*um3(x,y,z,t) +...
            2*pi*omega*sigm*(z + n1*uker1(x,y,z).*(x-y))/mum*cos(2*pi*omega*t);
        %u = ones(size(x));
    end
    function u = fp1(x,y,z,t)
        u = (-10*n2*uker2(x,y,z).*(z-y) +8*n2*(z.^3-y.^3+z.*y.^2+z.*x.^2-y.*x.^2-y.*z.^2) -...
             10*n2*uker1(x,y,z).*(z-y))*sin(2*pi*omega*t) - 4*pi^2*omega^2*epsp*up1(x,y,z,t) +...
             2*pi*omega*sigp*(x + n2*uker1(x,y,z).*uker2(x,y,z).*(y-z))/mup*cos(2*pi*omega*t);
        %u = ones(size(x));
    end
    function u = fp2(x,y,z,t)
        u = (-10*n2*uker2(x,y,z).*(x-z) +8*n2*(x.^3-z.^3+x.*z.^2+x.*y.^2-z.*y.^2-z.*x.^2) -...
             10*n2*uker1(x,y,z).*(x-z))*sin(2*pi*omega*t) - 4*pi^2*omega^2*epsp*up2(x,y,z,t) +...
             2*pi*omega*sigp*(y + n2*uker1(x,y,z).*uker2(x,y,z).*(z-x))/mup*cos(2*pi*omega*t);
        %u = ones(size(x));
    end
    function u = fp3(x,y,z,t)
        u = (-10*n2*uker2(x,y,z).*(y-x) +8*n2*(y.^3-x.^3+y.*x.^2+y.*z.^2-x.*z.^2-x.*y.^2) -...
             10*n2*uker1(x,y,z).*(y-x))*sin(2*pi*omega*t) - 4*pi^2*omega^2*epsp*up3(x,y,z,t) +...
             2*pi*omega*sigp*(z + n2*uker1(x,y,z).*uker2(x,y,z).*(x-y))/mup*cos(2*pi*omega*t);
        %u = ones(size(x));
    end

%% Diffusion coefficient function
    function u = Mu(x,y,z)
        u = Mum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Mup(x(id),y(id),z(id));
    end
    function u = Mum(x,y,z)
        u = mum*ones(size(x));
    end
    function u = Mup(x,y,z)
        u = mup*ones(size(x));
    end

%% Mass coefficient function
    function u = Epslon(x,y,z)
        u = Epslonm(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Epslonp(x(id),y(id),z(id));
    end
    function u = Epslonm(x,y,z)
        u = epsm*ones(size(x));
    end
    function u = Epslonp(x,y,z)
        u = epsp*ones(size(x));
    end
    
    function u = Sig(x,y,z)
        u = Sigm(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Sigp(x(id),y(id),z(id));
    end
    function u = Sigm(x,y,z)
        u = sigm*ones(size(x));
    end
    function u = Sigp(x,y,z)
        u = sigp*ones(size(x));
    end

%% Other function
    function u = one(x,y,z)
        u = ones(size(x));
    end
end