function pde = HodgeLaplacianFdata1
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.f = @f;
pde.exactu = @exactu;
pde.exactsigma = @exactsigma;
pde.gu = @exactu;
pde.gsigma = @exactsigma;
pde.gut = @gut;
pde.gdivu = 0;

    function s = f(p)
    x = p(:,1); y = p(:,2);
    s = [-sin(y), -cos(x)];
    end

    function s = exactu(p)
    x = p(:,1); y = p(:,2);
    s = [sin(y), cos(x)];
    end

    function s = exactsigma(p)
    x = p(:,1); y = p(:,2);
    s = -sin(x) - cos(y);   
    end

    function f = gut(p) % for unit square [0,1]^2
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2);
    u = exactu(p);
    leftbd = (abs(x)<eps);  % t = (0,-1); 
    f(leftbd) = - u(leftbd,2);
    rightbd = (abs(x-1)<eps); % t = (0,1); 
    f(rightbd) = u(rightbd,2);
    topbd = (abs(y-1)<eps);   % t = (-1,0)
    f(topbd) = - u(topbd,1);
    bottombd = (abs(y)<eps);% t = (1,0)
    f(bottombd) = u(bottombd,1);
    end

end