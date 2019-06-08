function pde = Stokesdata3
%% STOKESDATA3 data for Stokes equations
%
% A simple model of colliding flow. The force f = 0, the velocity  u1 =
% 1-y^2, u2 = 0, p = -2x.
%
% Parabolic inflow boundary condition, natural outflow boundary condition.
% The constant in p is chosen such that the Neumann boundary condition
% du_1/dx - p = on x = 1.
%
% Reference: page 233 in Finite Elements and Fast Iterative Solvers with
% Applications in Incompressible Fluid Dynamics. by Howard C. Elman, David
% J. Silvester, and Andrew J. Wathen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f', 0, 'exactp', @exactp, 'exactu', @exactu,'g_D',@g_D,'g_N',@g_N);

    % exact velocity
    function z = exactu(p)
    %x = p(:,1); 
    y = p(:,2);
    z(:,1) = 1-y.^2;
    z(:,2) = 0;
    end
    % exact pressure
    function z = exactp(p)
    x = p(:,1);% y = p(:,2);
    z = -2*x;
    end
    % Dirichlet boundary condition of velocity
    function z = g_D(p)
    z = exactu(p);
    end
    function z = g_N(p)
    z(:,1) = 2*ones(size(p,1),1);
    z(:,2) = zeros(size(p,1),1);
    end
end