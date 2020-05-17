function pde = arctandata
%% Arctan Circular Wave Front
% Reference: NIST AMR benchmark
% u = arctan(alpha*(r - r0))
% Neumann boundary left & bottom side
% Dirichlet boundary top & right side


pde = struct('f',@f,'d',@d,'exactu',@exactu,'g_D',@exactu,'g_N',@g_N,'Du',@Du);

x0 = -0.05;
y0 = -0.05;
r0 = 0.7;
alpha = 20;

    function z = f(p) % right hand side function
        x = p(:,1);
        y = p(:,2);
        r2 = (x - x0).^2 + (y - y0).^2;
        temp1 = alpha^3*(r2 - r0^2)-alpha;
        temp2 = alpha*(sqrt(r2) - r0);
        temp3 = sqrt(r2).*(temp2.^2 + 1).^2;
        z = temp1./temp3;
    end

    function z = d(p) % Diffusion constant
        z =  ones(size(p,1),1);
    end

    function z = exactu(p) % exact solution
        x = p(:,1);
        y = p(:,2);
        r2 = (x - x0).^2 + (y - y0).^2;
        z = atan(alpha*(sqrt(r2) - r0));
    end

    function z = Du(p)
        x = p(:,1);
        y = p(:,2);
        r2 = (x - x0).^2 + (y - y0).^2;
        temp = alpha*(sqrt(r2) - r0);
        z(:,1) = alpha*(x - x0)./sqrt(r2)./(temp.^2+1);
        z(:,2) = alpha*(y - y0)./sqrt(r2)./(temp.^2+1);
    end

    function s = g_N(p)
        Dup = Du(p);
        s = (p(:,1) == 1.0).*Dup(:,1) ... % on x=1, normal = (1,0)
            -(p(:,1) == -1.0).*Dup(:,1) + ... % on x = 0, norm = (-1,0)
            (p(:,2) == 1.0).*Dup(:,2) ... % on y=1, normal = (0,1)
            -(p(:,2) == -1.0).*Dup(:,2); % on y = 0, norm = (0,-1)
    end

end
