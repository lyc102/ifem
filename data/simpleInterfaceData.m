function pde = simpleInterfaceData
%%  SIMPLEINTERFACE data of the elliptic interface problem
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',0,'exactu',@exactu,'g_D',@exactu,'Du',@Du,'d',@d, 'phi', @phi);


r0 = pi/6.28;
phi = @(p) p(:,1).^2 + p(:,2).^2 - r0.^2;
betaMinus = 1;
betaPlus = 10;
alpha = 5;

    function K =  d(p)  % Diffusion constant
        idxPlus = (phi(p) >0);
        K = ones(size(p,1),1);
        K(idxPlus) = betaPlus;
    end

    function u = exactu(p)
        idxPlus = (phi(p) >0);
        idxMinus = 1 - idxPlus;
        r2 = p(:,1).^2 + p(:,2).^2;
        u = r2.^(alpha/2)./betaMinus.*idxMinus + ...
            (r2.^(alpha/2)./betaPlus + 1/betaMinus - 1/betaPlus).*idxPlus;
    end

    function Du = Du(p)
        idxPlus = (phi(p) >0);
        idxMinus = 1 - idxPlus;
        r2 = p(:,1).^2 + p(:,2).^2;
        Du(:,1) = 2*p(:,1).*r2.*(alpha/2-1)./betaMinus.*idxMinus + ...
            2*p(:,1).*r2.*(alpha/2-1)./betaPlus.*idxPlus;
        Du(:,2) = 2*p(:,2).*r2.*(alpha/2-1)./betaMinus.*idxMinus + ...
            2*p(:,2).*r2.*(alpha/2-1)./betaPlus.*idxPlus;
    end


end
