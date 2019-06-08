function pde = fracLapdata5
%% FRACLAPDATA5 f = D_x \chi_{x<0}
%
% f is a distribution: Dx of a Heaviside function of the half plane.
% 
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('intx',@intx); pde.g_N = 'intx';

    function z = intx(p)
        x1 = p(:,1); % x2 = p(:,2); % y = p(:,3);
        idx = (x1<0);
        z = zeros(size(p,1),1);
        z(idx) = 1;
    end

end