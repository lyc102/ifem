function pde = TorusTimeInitial4(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,a,b,intPt)
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
    'Sig',@Sig,'Sigm',@Sigm,'Sigp',@Sigp,...
    'E1',@E1,'E2',@E2,'E3',@E3,...
    'Et1',@Et1,'Et2',@Et2,'Et3',@Et3);

pde.mum = mum;
pde.mup = mup;
pde.sigm = sigm;
pde.sigp = sigp;
pde.epsm = epsm;
pde.epsp = epsp;
%% interface function
    function u = intf(x,y,z)
        u = (sqrt((x-x0).^2+(y-y0).^2)-r2).^2 + (z-z0).^2 - r1^2;
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
        u = zeros(size(x));
    end
    function u = um2(x,y,z,t)
        u = exp(-b*(a*(x-intPt)-omega*t).^2);
    end
    function u = um3(x,y,z,t)
        u = zeros(size(x));
    end
    function u = up1(x,y,z,t)
        u = zeros(size(x));
        %u = t*ones(size(x));
    end
    function u = up2(x,y,z,t)
        u = exp(-b*(a*(x-intPt)-omega*t).^2);
    end
    function u = up3(x,y,z,t)
        u = zeros(size(x));
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
        u = zeros(size(x));
    end
    function u = Dyum(x,y,z,t)
        u = zeros(size(x));
    end
    function u = Dzum(x,y,z,t)
        u = -exp(-b*(a*(x-intPt)-omega*t).^2)*2*a*b.*(a*(x-intPt)-omega*t);
    end
    function u = Dxup(x,y,z,t)
        u = zeros(size(x));
    end
    function u = Dyup(x,y,z,t)
        u = zeros(size(x));
    end
    function u = Dzup(x,y,z,t)
        u = -exp(-b*(a*(x-intPt)-omega*t).^2)*2*a*b.*(a*(x-intPt)-omega*t);
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
        u = zeros(size(x));
    end
    function u = fm2(x,y,z,t)
        u = zeros(size(x));
    end
    function u = fm3(x,y,z,t)
        u = zeros(size(x));
    end
    function u = fp1(x,y,z,t)
        u = zeros(size(x));
    end
    function u = fp2(x,y,z,t)
        u = zeros(size(x));
    end
    function u = fp3(x,y,z,t)
        u = zeros(size(x));
    end

%% Diffusion coefficient function
    function u = Mu(x,y,z)
        u = Mum(x,y,z);
        id = intf(x,y,z) > 0;
        u(id) = Mup(x(id),y(id),z(id));
    end
    function u = Mum(x,y,z)
        u = mum^(-1)*ones(size(x));
    end
    function u = Mup(x,y,z)
        u = mup^(-1)*ones(size(x));
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