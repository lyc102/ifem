function pde = fracLapdata9
%% FRACLAPDATA9 f = sin(2*pi*x1).*sin(pi*x2)
%
% sincos load
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f);

    function z = f(p)
        x1 = p(:,1); x2 = p(:,2); % y = p(:,3);
        z = sin(2*pi*x1).*sin(pi*x2);
    end

end