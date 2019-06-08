function pde = eddycurrentdata1
%% EDDYCURRENTDATA1 data for eddycurrent equation in two dimensions.
%
% See also: squareEddyCurrent
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.mu = 1;
pde.beta = 1;
pde.J = @J;
pde.exactu = @exactu;
pde.g_D = @exactu;
pde.curlu = [0 0];

    function s = J(p)
    s = exactu(p);
    end

    function s = exactu(p)
    x = p(:,1); y = p(:,2);
    s = [cos(x), -sin(y)];
    end
end