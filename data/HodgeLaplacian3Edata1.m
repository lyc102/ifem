function pde = HodgeLaplacian3Edata1
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.f = @f;
pde.u = @u;
pde.sigma = @sigma;
pde.gu = @u;
pde.gsigma = @sigma;
pde.gun = @gun;
pde.curlu = [0 0 0];
pde.gcurlu = 0;

    function s = f(p)
    x = p(:,1); y = p(:,2);
    s = [cos(x), -sin(y), zeros(size(p,1),1)];
    end

    function s = u(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = [cos(x), -sin(y), 2*z+1];
    end

    function s = sigma(p)
    x = p(:,1); y = p(:,2);
    s = sin(x)+cos(y)-2;   
    end

    function f = gun(p) % for unit cube [0,1]^3
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2); z = p(:,3);
    u = u(p);
    leftbd = (abs(x)<eps);    % n = (-1,0,0) 
    f(leftbd) = - u(leftbd,1);
    rightbd = (abs(x-1)<eps); % n = (1,0,0) 
    f(rightbd) = u(rightbd,1);
    frontbd = (abs(y)<eps);   % n = (0,-1,0)
    f(frontbd) = -u(frontbd,2);
    backbd = (abs(y-1)<eps);  % n = (0,1,0)
    f(backbd) = u(backbd,2);    
    topbd = (abs(z-1)<eps);   % n = (0,0,1)
    f(topbd) = u(topbd,3);
    bottombd = (abs(z)<eps);  % n = (0,0,-1)
    f(bottombd) = - u(bottombd,3);
    end


end