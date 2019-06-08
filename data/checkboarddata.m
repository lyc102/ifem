function pde = checkboarddata
%% CHECKBOARDPDATA f = (x1.^2 - 1).*(x2.^2 - 1), d = 161.4476387975881
%
% Data for jump coefficients (checkboard) example.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'d',@d);

    function z = f(p)
        x1 = p(:,1); x2 = p(:,2); % y = p(:,3);
        z = (x1.^2 - 1).*(x2.^2 - 1);
    end

    function K =  d(p)  % Diffusion constant
        idx = (p(:,1).*p(:,2) >0);
        K  = ones(size(p,1),1);
        K(idx) = 161.4476387975881;
    end
end