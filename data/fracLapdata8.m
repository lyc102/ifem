function pde = fracLapdata8
%% FRACLAPDATA8 f = D_x \chi_{|x|_{\infty}<0.5}
%
% f is a distribution Dx of characteristic function of a square.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('intx',@intx); pde.f = 'intx';

    function z = intx(p)
        idx = (max(abs(p),[],2) < 0.5);
        z = zeros(size(p,1),1);        
        z(idx) = 1;
    end

end