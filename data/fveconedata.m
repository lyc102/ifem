function pde = fveconedata
%% FVECONEDATA f = [1 1] or [1 1 1]
% 
% f = [1 1] or [1 1 1].
%
% No exact solution.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f);

    function z = f(p) %#ok<*INUSD>
        z = ones(size(p,1),size(p,2));
    end

end