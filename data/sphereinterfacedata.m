function pde = sphereinterfacedata(r, betaPlus, betaMinus)

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du,'phi',@phi, 'exactw',@exactw,...
    'exactq',@exactq, 'd',@d,'exactuplus',@exactuplus,'exactuminus',@exactuminus,...
    'Duplus',@Duplus,'Duminus',@Duminus,'dplus',@dplus,'dminus',@dminus);

% load data (right hand side function)

    function s = dplus(p)
        s = betaPlus;
    end
    function s = dminus(p)
        s = betaMinus;
    end
    function s = d(p)
        N = size(p,1);
        s = zeros(N,1);
        isInside = msign(phi(p)) <= 0;
        s(isInside) = betaMinus;
        s(~isInside) = betaPlus;
    end
    function s = f(p)
        isInside = msign(phi(p))<=0;
        s = zeros(size(p,1),1);
        s(isInside) = - 6*betaMinus*sum(p(isInside,:),2);
        s(~isInside) = - 6*betaPlus*sum(p(~isInside,:),2);
    end
% exact solution
    function s = exactu(p)
        N = size(p,1);
        s = zeros(N,1);
        isInside = msign(phi(p)) <= 0;
        s(isInside) = exactuminus(p(isInside,:));
        s(~isInside) = exactuplus(p(~isInside,:));
    end

    function s = exactuplus(p)
        s = sum(p.^3,2);
    end
    function s = exactuminus(p)
        s = sum(p.^3,2);
    end

% Dirichlet boundary condition
    function s = g_D(p)
        s = exactuplus(p);
    end
% Derivative of the exact solution
    function s = Du(p)
        N = size(p,1);
        s = zeros(N,3);
        isInside = msign(phi(p)) <= 0;
        s(isInside,:) = Duminus(p(isInside,:));
        s(~isInside,:) = Duplus(p(~isInside,:));
    end

    function s = Duminus(p)
        s = zeros(size(p));
        s(:,1) = 3*p(:,1).^2;
        s(:,2) = 3*p(:,2).^2;
        s(:,3) = 3*p(:,3).^2;
    end

    function s = Duplus(p)
        s = zeros(size(p));
        s(:,1) = 3*p(:,1).^2;
        s(:,2) = 3*p(:,2).^2;
        s(:,3) = 3*p(:,3).^2;
    end

    function s = exactw(p)
        s = exactuplus(p) - exactuminus(p);
    end

    function s = exactq(p)
        % jump of the flux [\beta du/dn]
         x = p(:,1); y = p(:,2); z =p(:,3);
        l = sqrt(sum(p.^2,2));
        n = [x./l,y./l, z./l];
        p = r*n;
        z1 = Duplus(p);
        z2 = Duminus(p);
        s = [betaPlus*z1(:,1)-betaMinus*z2(:,1), betaPlus*z1(:,2)-betaMinus*z2(:,2),...
             betaPlus*z1(:,3)-betaMinus*z2(:,3)];
        s = sum(s.*n,2);
        
    end

    function s = phi(p) % level set function
        x = p(:,1); y = p(:,2); z =p(:,3);
        s = x.^2 + y.^2 +z.^2 - r^2;
    end

end