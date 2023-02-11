function pde = TorusCenterTime(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,r,stre)
%% USAGE: polynomial solution for Poisson equation
%  Last Modified: 02/21/2020 by Xu Zhang

%% PDE Structure
pde = struct('intf',@intf,...
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
        u = 2*pi*omega*sin(2*pi*omega*t).*exp(-L2Norm(x,y,z)/r^2).*x*stre;
    end
    function u = fm2(x,y,z,t)
        u = 2*pi*omega*sin(2*pi*omega*t).*exp(-L2Norm(x,y,z)/r^2).*y*stre;
    end
    function u = fm3(x,y,z,t)
        u = 2*pi*omega*sin(2*pi*omega*t).*exp(-L2Norm(x,y,z)/r^2).*z*stre;
    end
    function u = fp1(x,y,z,t)
        u = 2*pi*omega*sin(2*pi*omega*t).*exp(-L2Norm(x,y,z)/r^2).*x*stre;
    end
    function u = fp2(x,y,z,t)
        u = 2*pi*omega*sin(2*pi*omega*t).*exp(-L2Norm(x,y,z)/r^2).*y*stre;
    end
    function u = fp3(x,y,z,t)
        u = 2*pi*omega*sin(2*pi*omega*t).*exp(-L2Norm(x,y,z)/r^2).*z*stre;
    end
    function r = L2Norm(x,y,z)        
        r = x.^2+y.^2+z.^2;       
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