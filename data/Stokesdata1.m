function pde = Stokesdata1
%% STOKESDATA1 data for Stokes equations
%
% A simple model of colliding flow. The force f = 0, the velocity  u1 =
% 20xy^3, u2 = 5x^4-5y^4, p = 60x^2y - 20y^3.
%
% Dirichlet boundary condition is imposed.
%
% Reference: page 237 in Finite Elements and Fast Iterative Solvers with
% Applications in Incompressible Fluid Dynamics. by Howard C. Elman, David
% J. Silvester, and Andrew J. Wathen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f', 0, 'exactp', @exactp, 'exactu', @exactu,'g_D',@exactu, 'exactw', @exactw);

    % exact velocity
    function z = exactu(p)
        x = p(:,1); y = p(:,2);
        z(:,1) = 20*x.*y.^3;
        z(:,2) = 5*x.^4-5*y.^4;
    end
    % exact pressure
    function z = exactp(p)
        x = p(:,1); y = p(:,2);
        z = 60*x.^2.*y-20*y.^3-5;
    end
    % exact vorticity
    function z = exactw(p)
        x = p(:,1); y = p(:,2);
        z = 20*x.^3-60*x.*y.^2;
    end
end
