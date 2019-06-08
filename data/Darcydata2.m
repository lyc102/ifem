function pde = Darcydata2
%% Darcydata2 anisotropic tensor
%
% p = x^3y+y^4 + sin(pi*x)cos(pi*y)
% K = [1 + 4*(x^2 + y^2)    3*x*y]
%   = [3*x*y   1 + 11*(x^2 + y^2)]

pde = struct('exactp',@exactp,'K',@K,'exactu',@exactu,'Dp',@gradp,...
             'f', @f,'g_D',@g_D,'g_N',@g_N);

    function s = K(pt)
       x = pt(:,1); y = pt(:,2);
%        s(:,3) = 3*x.*y;        
       s(:,1) = 1 + 4*(x.^2 + y.^2);
       s(:,2) = 1 + 11*(x.^2 + y.^2);
    end
    function s = exactp(pt)
       x = pt(:,1); y = pt(:,2);
%        s = x.^3.*y + y.^4 + sin(pi*x).*cos(pi*y);
       s = x.^2 + y.^2;
    end
    function s = gradp(pt)
       x = pt(:,1); y = pt(:,2);
%        s(:,2) = x.^3 + 4*y.^3 - pi*sin(pi*x).*sin(pi*y);
%        s(:,1) = 3*x.^2.*y + pi*cos(pi*x).*cos(pi*y);       
       s(:,2) = 2*y;
       s(:,1) = 2*x;       
    end
    function s = exactu(pt)
       Dp = gradp(pt);  % u = K*grad(p)
       K = K(pt);
       s(:,2) = K(:,3).*Dp(:,1) + K(:,2).*Dp(:,2);
       s(:,1) = K(:,1).*Dp(:,1) + K(:,3).*Dp(:,2);
    end
    function s = f(pt)
        x = pt(:,1); y = pt(:,2);
%         s = 2+24*x.^2+14*y.^2+6*x.^2+2+22*x.^2+66*y.^2;
        s = 2+24*x.^2+16*y.^2++2+22*x.^2+66*y.^2;
        s = -s;
    end
    function s = g_D(pt)
       s = exactp(pt); 
    end
    function f = g_N(p,vargin)
        if nargin > 1
            f = dot(exactu(p),vargin,2);
        else
            f = zeros(size(p,1),1);
            x = p(:,1); 
            y = p(:,2);
            u = exactu(p);
            leftbd = (abs(x)<eps);  % n = (-1,0); 
            f(leftbd) = - u(leftbd,1);
            rightbd = (abs(x-1)<eps); % n = (1,0); 
            f(rightbd) = u(rightbd,1);
            topbd = (abs(y-1)<eps);   % n = (0,1)
            f(topbd) = u(topbd,2);
            bottombd = (abs(y)<eps);% n = (0,-1)
            f(bottombd) = - u(bottombd,2);    
        end
    end
end