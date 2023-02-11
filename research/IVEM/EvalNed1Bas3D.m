function f = EvalNed1Bas3D(bas_total, tID, x, y, z, bas_ind, dind, vind)

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
    
    if vind == 1
        bas = bas_total.bas_curl_x(tID,:,bas_ind);
    elseif vind == 2
        bas = bas_total.bas_curl_y(tID,:,bas_ind);
    elseif vind == 3
        bas = bas_total.bas_curl_z(tID,:,bas_ind);
    end
    
    c1 = bas(:,1); c2 = bas(:,2); c3 = bas(:,3); c4 = bas(:,4);
    
    f = c1.*one + c2.*x + c3.*y + c4.*z;
    
elseif dind == 1
    
    bas_x = bas_total.bas_curl_x(tID,:,bas_ind);
    bas_y = bas_total.bas_curl_y(tID,:,bas_ind);
    bas_z = bas_total.bas_curl_z(tID,:,bas_ind);
    
    c12 = bas_x(:,2); c13 = bas_x(:,3); c14 = bas_x(:,4);
    c22 = bas_y(:,2); c23 = bas_y(:,3); c24 = bas_y(:,4);
    c32 = bas_z(:,2); c33 = bas_z(:,3); c34 = bas_z(:,4);
    
    if vind == 1
        f = c33 - c24;
    elseif vind == 2
        f = c14 - c32;
    elseif vind == 3
        f = c22 - c13;
    end
    
end

end