function pde = StokesZulehnerdata
%% STOKESDATA2 data for Stokes equations
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f', @f, 'exactp', @exactp, 'exactu', @exactu,'g_D',@g_D);

    function z = f(p)
        x = p(:,1); y = p(:,2);
        z(:,1) = zeros(size(x));
        z(:,2) = 4*cos(x).*cos(y);
    end
    % exact velocity
    function z = exactu(p)
        x = p(:,1); y = p(:,2);
        z(:,1) = sin(x).*sin(y);
        z(:,2) = cos(x).*cos(y);
    end
    % exact pressure
    function z = exactp(p)
        x = p(:,1); y = p(:,2);
        z = 2*cos(x).*sin(y) - 2*sin(1)*(1-cos(1));
    end
    % Dirichlet boundary condition of velocity
    function z = g_D(p)
        z = exactu(p);
    end
end