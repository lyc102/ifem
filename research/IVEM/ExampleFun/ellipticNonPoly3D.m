function pde = ellipticNonPoly3D(aa,bb,cc)
%% USAGE: Non-polynomial solution for Poisson equation
%
%     u = (x^2+y^2+z^2)^(5/2);
%     f = -30*(x^2+y^2+z^2)^(3/2);
%     Du = (5*x.*(x.^2+y.^2).^(3/2); 5*y.*(x.^2+y.^2).^(3/2);
%           5*z.*(x.^2+y.^2).^(3/2));
%
% Last Modified: 05/24/2016 by Xu Zhang
% Last Modified: 07/02/2020 by Xu Zhang

%%
pde = struct('f',@f,'exactu',@exactu,'gD',@gD,'Du',@Du,'A',@A,'b',@b,'b1',@b1,...
    'b2',@b2,'b3',@b3,'c',@c,'one',@one,'Dxu',@Dxu,'Dyu',@Dyu,'Dzu',@Dzu);

% exact solution
    function u = exactu(x,y,z)
        u = (x.^2+y.^2+z.^2).^(5/2);
    end
% right hand side function
    function u = f(x,y,z)
        u = -30*A(x,y,z).*(x.^2+y.^2+z.^2).^(3/2) + ...
            b1(x,y,z).*Dxu(x,y,z) + b2(x,y,z).*Dyu(x,y,z) + ...
            b3(x,y,z).*Dzu(x,y,z) + c(x,y,z).*exactu(x,y,z);
    end
% Dirichlet boundary condition
    function u = gD(x,y,z)
        u = exactu(x,y,z);
    end
% Derivative of the exact solution
    function u = Dxu(x,y,z)
        u = 5*x.*(x.^2+y.^2+z.^2).^(3/2);
    end
    function u = Dyu(x,y,z)
        u = 5*y.*(x.^2+y.^2+z.^2).^(3/2);
    end
    function u = Dzu(x,y,z)
        u = 5*z.*(x.^2+y.^2+z.^2).^(3/2);
    end
    function u = Du(x,y,z)
        ux = Dxu(x,y,z);
        uy = Dyu(x,y,z);
        uz = Dzu(x,y,z);
        u = [ux;uy;uz];
    end
% Diffusion coefficient function
    function u = A(x,y,z)
        u = aa*ones(size(x)).*ones(size(y)).*ones(size(z));
    end
% Advection coefficient function
    function u = b(x,y,z)
        ux = b1(x,y,z);
        uy = b2(x,y,z);
        uz = b3(x,y,z);
        u = [ux;uy;uz];
    end
    function u = b1(x,y,z)
        u = bb(1)*ones(size(x)).*ones(size(y)).*ones(size(z));
    end
    function u = b2(x,y,z)
        u = bb(2)*ones(size(x)).*ones(size(y)).*ones(size(z));
    end
    function u = b3(x,y,z)
        u = bb(3)*ones(size(x)).*ones(size(y)).*ones(size(z));
    end
% Reaction coefficient function
    function u = c(x,y,z)
        u = cc*ones(size(x)).*ones(size(y)).*ones(size(z));
    end
% one function
    function u = one(x,y,z)
        u = ones(size(x)).*ones(size(y)).*ones(size(z));
    end
end