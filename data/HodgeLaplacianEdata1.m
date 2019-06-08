function pde = HodgeLaplacianEdata1
%% HODGELAPLACIANEDATA1
%
% u = [cos(x), -sin(y)];
% sigma = -div u = sin(x) + cos(y);
% f = -grad div u + curl rot u = - Lap u;
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.f = @f;
pde.u = @exactu;
pde.sigma = @exactsigma;
pde.gu = @exactu;
pde.gsigma = @exactsigma;
pde.gun = @gun;
pde.grotu = 0;

    function s = f(p)
    s = exactu(p);
    end

    function u = exactu(p)
    x = p(:,1); y = p(:,2);
    u = [cos(x), -sin(y)];
    end

    function s = exactsigma(p)
    x = p(:,1); y = p(:,2);
    s = sin(x) + cos(y);   
    end

    function s = gun(p) % for unit square [0,1]^2
    s = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2);
    u = exactu(p);
    leftbd = (abs(x)<eps);  % n = (-1,0); 
    s(leftbd) = - u(leftbd,1);
    rightbd = (abs(x-1)<eps); % n = (1,0); 
    s(rightbd) = u(rightbd,1);
    topbd = (abs(y-1)<eps);   % n = (0,1)
    s(topbd) = u(topbd,2);
    bottombd = (abs(y)<eps);% n = (0,-1)
    s(bottombd) = - u(bottombd,2);
    end

end