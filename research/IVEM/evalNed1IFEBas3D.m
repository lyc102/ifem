function f = evalNed1IFEBas3D(bas_total, x, y, z, dind, v_ind)

%% USAGE: evaluate P1 solution or interpolation at certain point(s).
%
% INPUTS:
% bas --- coefficient of basis [c1,c2,c3,c4], ci can be vectors
%            c1 + c2*x + c3*y + c4*z
% x, y, z --- coordinates of the points to evaluate can becolumn vectors,
% dind --- derivative info for FE function
%            d = [0,0,0]: function value
%            d = [1,0,0]: Dx value
%            d = [0,1,0]: Dy value
%            d = [0,0,1]: Dz value
% OUTPUTS:
% f --- the value of the FE/DG solution or interpolation at point(x,y).
%
% Last Modified: 06/19/2020 by Xu Zhang
% Last Modified: 07/02/2020 by Xu Zhang
%%
one = ones(size(x));

if dind == 0
     
    a1 = bas_total(:,1); a2 = bas_total(:,2); a3 = bas_total(:,3); 
    b1 = bas_total(:,4); b2 = bas_total(:,5); b3 = bas_total(:,6);
    
    if v_ind == 1
        f = b1.*one - a3.*y + a2.*z;
    elseif v_ind == 2
        f = b2.*one + a3.*x - a1.*z;
    elseif v_ind == 3
        f = b3.*one - a2.*x + a1.*y;
    end
    
elseif dind == 1
    
    a1 = bas_total(:,1); a2 = bas_total(:,2); a3 = bas_total(:,3); 
    
    if v_ind == 1
        f = 2*a1*ones(1,size(x,2));
    elseif v_ind == 2
        f = 2*a2*ones(1,size(x,2));
    elseif v_ind == 3
        f = 2*a3*ones(1,size(x,2));
    end
    
end

end