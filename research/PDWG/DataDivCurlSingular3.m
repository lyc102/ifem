function pde = DataDivCurlSingular3
%%
%           curl u = g   in \Omega,
%           div (eps* u) = f    in \Omega,
%           n \cdot (eps * u) = g_N  on \Gamma_0
% Reference: modified based on Example 2 in
% A weak Galerkin least-squares finite element method for div-curl systems

pde.Eps = @(p) [1+p(:,1).*0 0.*p(:,1) 0.*p(:,1); ...
    0*p(:,2) 1+p(:,2).*0 0*p(:,2); ...
    0*p(:,2) 0*p(:,2) 1+p(:,3).*0];
pde.exactu = @exactu;
pde.g = @curlu;
pde.f = @divu;
pde.curlu = @curlu;

    function s = exactu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        theta = atan2(y,x);
%         theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        r2 = x.^2+y.^2;
        s(:,1) = x.*(1-x);
        s(:,2) = y.*(1-y);
        s(:,3) = r2.^(1/3).*sin(2*theta).*z.*(1-z);
    end

    function s = curlu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        theta = atan2(y,x);
        r2 = x.^2+y.^2;

        s(:,1) = z.*(1-z).*(2*y.*sin(2*theta) + 6*x.*cos(2*theta))./r2.^(2/3)/3;
        s(:,2) = z.*(1-z).*(-2*x.*sin(2*theta) + 6*y.*cos(2*theta))./r2.^(2/3)/3;
        s(:,3) = z*0;
    end

    function s = divu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        theta = atan2(y,x);
        r2 = x.^2+y.^2;
        s = 2 - 2*x - 2*y + (1-2*z).*r2.^(1/3).*sin(2*theta);
    end


end