function pde = Stokesdata2
%% STOKESDATA2 data for Stokes equations
%
% The solution u is a polynomial satisfying the zero Dirichlet boundary
% condition. 
%
% Created by Ming Wang.
%
pde = struct('f', @f, 'exactp', @exactp, 'exactu',@exactu, 'exactw', @exactw, 'Dw',@exactDw,'g_D',@exactu);
%
%% subfunction
    function z = f(p) % load data (right hand side function)
        x = p(:,1); y = p(:,2);
        z(:,1) = 2^10*((1-6*x+6*x.^2).*(y-3*y.^2+2*y.^3) ...
            + (x.^2-2*x.^3+x.^4).*(-3+6*y)...
            - (-3+6*x).*(y.^2-2*y.^3+y.^4) );
        z(:,2) = -(2^10)*( (-3+6*x).*(y.^2-2*y.^3 ...
            + y.^4)+(x-3*x.^2+2*x.^3).*(1-6*y+6*y.^2)...
            + (1-6*x+6*x.^2).*(y-3*y.^2+2*y.^3) );
    end

    function z = exactu(p)
        x = p(:,1); y = p(:,2);
        z(:,1) = -(2^8)*(x.^2-2*x.^3+x.^4).*(2*y-6*y.^2+4*y.^3);
        z(:,2) = 2^8*(2*x-6*x.^2+4*x.^3).*(y.^2-2*y.^3+y.^4);
    end

    function z = exactp(p)
        x = p(:,1); y = p(:,2);
        z = -(2^8)*(2-12*x+12*x.^2).*(y.^2-2*y.^3+y.^4);
    end

    function z = exactw(p)
        x = p(:,1); y = p(:,2);
        z = (y.^4 - 2*y.^3 + y.^2).*(3072*x.^2 - 3072*x + 512) + ...
            (256*x.^4 - 512*x.^3 + 256*x.^2).*(12*y.^2 - 12*y + 2);
    end

    function z = exactDw(p)
        x = p(:,1); y = p(:,2);
        z(:,1) = (y.^4 - 2*y.^3 + y.^2).*(2*3072*x - 3072) + ...
            (4*256*x.^3 - 3*512*x.^2 + 2*256*x).*(12*y.^2 - 12*y + 2);
        z(:,2) = (4*y.^3 - 6*y.^2 + 2*y).*(3072*x.^2 - 3072*x + 512) + ...
                 (256*x.^4 - 512*x.^3 + 256*x.^2).*(24*y - 12);
    end
end