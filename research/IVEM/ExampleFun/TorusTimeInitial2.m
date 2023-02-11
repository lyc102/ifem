function pde = TorusTimeInitial2(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,a,b,intPt)
%% USAGE: polynomial solution for Poisson equation
%  Last Modified: 02/21/2020 by Xu Zhang

%% PDE Structure
pde = struct('intf',@intf,...
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

%%%%% intial condition
    function u = E1(x,y,z,t)
        u = exp(-b*(a*(x-intPt)-omega*t).^2);
%         id = pulregion(x,y,z) > 0;
%         u(id) = zeros(size(x(id)));
    end
    function u = E2(x,y,z,t)
        u = zeros(size(x));
    end
    function u = E3(x,y,z,t)
        u = zeros(size(x));
    end
    function u = Et1(x,y,z,t)
        u = exp(-b*(a*(x-intPt)-omega*t).^2)*2*omega*b.*(a*(x-intPt)-omega*t);
    end
    function u = Et2(x,y,z,t)
        u = zeros(size(x));
    end
    function u = Et3(x,y,z,t)
        u = zeros(size(x));
    end
%     function u = pulregion(x,y,z)
%        u = (x-intPt).^2/(0.15)^2 - 1; 
%     end
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