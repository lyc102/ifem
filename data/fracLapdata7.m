function pde = fracLapdata7
%% FRACLAPDATA7 f = \chi_{|x|_{\infty}<0.5}
%
% Piecewise constant load: characterisitic function of a square.

% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f);

    function z = f(p)
        idx = (max(abs(p),[],2) < 0.5);
        z = zeros(size(p,1),1);        
        z(idx) = 1;
    end

end