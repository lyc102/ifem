function pde = tensordata1
%% TENSORDATA1 data for interface problem
%
% Created by Huayi Wei Monday, 27 June 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du,'d',@Diff,'dinv',@dinv);

% load data (right hand side function)
% [10^4 0; 0 1] (0<x<1/2), [1 0; 0 2] (1/2<x<1);
%lmat = [10^4 0; 0 1]; rmat = [1 0; 0 2]; % matrix on both sides of x = 1/2;
lmat = 2*[1 0; 0 1]; rmat = 2*[1 0; 0 1]; % matrix on both sides of x = 1/2;
    function z = f(p)
    x = p(:,1); y = p(:,2);
    z = -(2*lmat(1,1)*(y.^2-y)+2*lmat(2,2)*(x.^2-x)).*(x<1/2) ...
        -(2*rmat(1,1)*(y.^2-y)+2*rmat(2,2)*(x.^2-x)).*(x>1/2);   
    end
    %exact solution
    function z = exactu(p)
    x = p(:,1); y = p(:,2);
    z = (x.^2-x).*(y.^2-y);
    end
    % Dirichlet boundary condition
    function z = g_D(p)
    z = exactu(p);
    end
    % Derivative of the exact solution
    function z = Du(p) % maybe need to change, let's see.
    x = p(:,1); y = p(:,2);
    z(:,1) = (2*x-1).*(y.^2-y);
    z(:,2) = (2*y-1).*(x.^2-x);
    end
    function z = Diff(p)
    x = p(:,1);
    NT = length(x);
    z = zeros(NT,2,2);
    z(:,1,1) = (x<1/2)*lmat(1,1) + (x>1/2)*lmat(2,2);
    z(:,2,2) = (x<1/2)*rmat(1,1) + (x>1/2)*rmat(2,2);
    end

    function z = dinv(p)
    x = p(:,1);
    NT = length(x);
    lmat = lmat^(-1); rmat = rmat^(-1);
    z = zeros(NT,2,2);
    z(:,1,1) = (x<1/2)*lmat(1,1) + (x>1/2)*lmat(2,2);
    z(:,2,2) = (x<1/2)*rmat(1,1) + (x>1/2)*rmat(2,2);
    end
end
