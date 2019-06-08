function pde = fracLapdata6
%% FRACLAPDATA6 f = D_x \chi_{B_1/2} 
%
% f is a distribution Dx of characteristic function of ball with radius 1/2.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('intx',@intx); pde.f = 'intx';

    function z = intx(p)
        x1 = p(:,1);  x2 = p(:,2); % y = p(:,3);
        idx = ((x1.^2 + x2.^2) < 0.25);
        z = zeros(size(p,1),1);
        z(idx) = 1;
    end

end