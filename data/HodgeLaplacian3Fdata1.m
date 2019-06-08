function pde = HodgeLaplacian3Fdata1
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.f = @f;
pde.exactu = @exactu;
pde.exactsigma = 0;
pde.gu = @exactu;
pde.gsigma = 0;
pde.gu = @exactu;
pde.curlu = 0;
pde.gdivu = @gdivu;

    function s = f(p)
    x = p(:,1); y = p(:,2);
    s = [cos(x), -sin(y), zeros(size(p,1),1)];
    end

    function s = exactu(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = [cos(x), -sin(y), 2*z+1];
    end

    function s = gdivu(p)
    x = p(:,1); y = p(:,2);
    s = sin(x)+cos(y)-2;
    end


end