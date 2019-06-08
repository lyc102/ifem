function pde = jumpmgdata1
%% Data of JUMPMG1
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'g_D',@g_D,'d',@omega);

    function s = f(p) % load data (right hand side function)
    s = zeros(size(p,1),1);
    end

    function s = g_D(p) % Dirichlet boundary condition
    s = zeros(size(p,1),1);
    x = p(:,1); 
    idx = (abs(x-1)<eps);
    s(idx) = 1;
    s(~idx) = 0;
    end

    function c = omega(p) % diffusion constant
    global epsilon
    c = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2); z = p(:,3);
    idx = (x>0) & (x<1) & (y>0) & (y<1) & (z>0) & (z<1);
    c(idx) = 1/epsilon;
    c(~idx) = 1;
    end
end