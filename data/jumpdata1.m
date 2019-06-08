function pde = jumpdata1
%% JUMPDATA1 data for interface problem
%
% Created by Huayi Wei Monday, 27 June 2011.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du,'d',@Diff);

% load data (right hand side function)
ljump = 1000; rjump = 1; 
    function z = f(p)
    x = p(:,1); y = p(:,2);
    z = -2*ljump*(y.^2-y+x.^2-x).*(x<1/2) ...
        -2*rjump*(y.^2-y+x.^2-x).*(x>1/2);
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
    %NT = length(x);
    %z = zeros(NT,1);
    z=ljump*(x<1/2) + rjump*(x>1/2);
    end
end
