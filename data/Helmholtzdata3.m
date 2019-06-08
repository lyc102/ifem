function pde = Helmholtzdata3
%% HELMHOTLZDATA4 trigonometric  data for PML Helmholtz equation
% -laplace u -k^2 u = f;
% The Intresting domain is [0,1]^2; The computing domain including the PML
%domain is [-0.1,1.1]^2;
% boundary condition
% u =0 on the PML boundary; 
% Created by Jie Zhou.
% The data come from
% Advances in Iterative Methods and Preconditioners for the Helmholtz Equation
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'k2',@k2,'g_D',@g_D,'d',@d);

    % load data (right hand side function)
    function rhs =  f(p)

      x = p(:,1); y = p(:,2);
     k0 = size(p,1);
    rhs = zeros(k0,1);
     ii = find(abs(x-0.5)<=1.0e-2&abs(y-0.5)<=1.0e-2);
    rhs(ii) = 10000;
    end



    function PMLp = d(p)
    global omegal k  c
    c = 340;
    omegal = k*340;
          i = sqrt(-1);
          x = p(:,1); y = p(:,2);
      size0 = size(x,1);
     sigma1 = ones(size0,1);
     sigma2 = ones(size0,1);
    index = find((x<0)|(x>1)|(y>1)|(y<0)); 
    sigma1(index) = x(index);
    sigma2(index) = y(index);
    
    s1 = 1 + sigma1/((i)*omegal);
    s2 = 1 + sigma2/((i)*omegal);
    
    PMLp = [s2./s1 s1./s2];
    end
    function s =  g_D(p)
    s = zeros(size(p,1),1);
    end
    function wavenumber = k2(p)
    global k omegal  c
    omegal = k*340;
          i = sqrt(-1);
          x = p(:,1); y = p(:,2);
      size0 = size(x,1);
     sigma1 = ones(size0,1);
     sigma2 = ones(size0,1);
    index = find((x<0)|(x>1)|(y>1)|(y<0)); 
    sigma1(index) = x(index);
    sigma2(index) = y(index);
    
    s1 = 1 + sigma1/((i)*omegal);
    s2 = 1 + sigma2/((i)*omegal);
    wavenumber = k.^2.*s1.*s2;
    end





end
