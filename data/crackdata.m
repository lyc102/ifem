function pde = crackdata
%% CRACKDATA data of Crack Problem
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@exactu,'Du',@Du);


    function z = f(p)   % load data (right hand side function)
    z = ones(size(p,1),1);
    end

    function u = exactu(p) % exact solution
    r = sqrt(sum(p.^2,2));
    u = sqrt(0.5*(r-p(:,1)))-0.25*r.^2;
    end

    function uprime = Du(p) % derivative of the exact solution
    r = sqrt(sum(p.^2,2));
    uprime(:,1) = (p(:,1)./r-1)./sqrt(8*(r-p(:,1)))-0.5*p(:,1); % u_x
    uprime(:,2) = p(:,2)./r./sqrt(8*(r-p(:,1)))-0.5*p(:,2);     % u_y
    end
end