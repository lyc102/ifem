function pde = cubeafemdata
%% CUBEAFEMDATA Data of an example of AFEM in a cube
%
%   f = 20*exp(-10*(r^2)).*(3-20*r^2);
%   u = exp(-10*r^2);
%  Du = -20e^(-10*r^2)(x,y,z); 
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du);

    function s = f(p) % load data (right hand side function)
    r2 = sum(p.^2,2);
    s = 20*exp(-10*(r2)).*(3-20*r2);
    end

    function s = g_D(p) % Dirichlet boundary condition
    s = exactu(p);
    end

    function s = exactu(p) % exact solution
    r2 = sum(p.^2,2);
    s = exp(-10*r2);
    end

    function s = Du(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    r2 = sum(p.^2,2);
    s(:,1) = -20*exp(-10*r2).*x;
    s(:,2) = -20*exp(-10*r2).*y;
    s(:,3) = -20*exp(-10*r2).*z;
    end
end