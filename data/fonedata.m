function pde = fonedata
%% FONEDATA f = 1
% 
% f = 1.
%
% No exact solution.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f);

    function z = f(p) %#ok<*INUSD>
        z = 1;
    end

end