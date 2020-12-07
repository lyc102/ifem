function pde = simpleInterfaceData(varargin)
%%  SIMPLEINTERFACE data of the elliptic interface problem
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@exactu,'Du',@Du,'d',@d, 'phi', @phi);

if isempty(varargin)
    betaMinus = 1;
    betaPlus = 5;
else
    betaMinus = varargin{1};
    betaPlus = varargin{2};
end
r0 = pi/6;
phi = @(p) p(:,1).^2 + p(:,2).^2 - r0.^2;

alpha = 5;

    function K =  d(p)  % Diffusion constant
        idxPlus = (phi(p) >0);
        K = betaMinus*ones(size(p,1),1);
        K(idxPlus) = betaPlus;
    end

    function u = exactu(p)
        idxPlus = (phi(p) >0);
        idxMinus = 1 - idxPlus;
        r2 = p(:,1).^2 + p(:,2).^2;
        u = r2.^(alpha/2)./betaMinus.*idxMinus + ...
            (r2.^(alpha/2)./betaPlus + r0^alpha*(1/betaMinus - 1/betaPlus)).*idxPlus;
    end

    function uprime = Du(p)
        idxPlus = (phi(p) >0);
        idxMinus = 1 - idxPlus;
        r2 = p(:,1).^2 + p(:,2).^2;
        uprime(:,1) = alpha*(p(:,1).*r2.*(alpha/2-1)./betaMinus.*idxMinus + ...
            p(:,1).*r2.*(alpha/2-1)./betaPlus.*idxPlus);
        uprime(:,2) = alpha*(p(:,2).*r2.*(alpha/2-1)./betaMinus.*idxMinus + ...
            p(:,2).*r2.*(alpha/2-1)./betaPlus.*idxPlus);
    end

    function s = f(p)
        r2 = p(:,1).^2 + p(:,2).^2;
        s = -alpha^2*r2.^(alpha/2-1);
    end


end
