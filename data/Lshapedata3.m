function pde = Lshapedata3
%% LSHAPEDATA data of Lshape Problem
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


pde = struct('f',0,'exactu',@exactu,'g_D',@exactu,'Du',@Du);

        
    function u = exactu(p) % exact solution
    x = p(:,1); y = p(:,2);
    r = sqrt(x.^2+y.^2);
    theta = atan2(p(:,2),p(:,1));
    theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
    u = r.^(2/3).*sin(2*theta/3);
    end

    function uprime = Du(p) % exact solution
    x = p(:,1); y = p(:,2);
    r = sqrt(x.^2+y.^2);
    theta = atan2(y,x);
    theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
    uprime(:,1) = 2/3*r.^(-1/3).*sin(2*theta/3).*x./r ...
                - 2/3*r.^(2/3).*cos(2*theta/3).*y./r.^2;
    uprime(:,2) = 2/3*r.^(-1/3).*sin(2*theta/3).*y./r ...
                + 2/3*r.^(2/3).*cos(2*theta/3).*x./r.^2;
    uprime(:,3) = 0;
    end
end