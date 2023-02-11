function gw = gaussWtetra(ng)

%% USAGE: generate Gaussian weights on a tetrahedron
%
% INPUTS:
% ngauss --- number of Gaussian nodes in a tetrahedron
%            possible values = 1, 4, 5, 11, 15.
%
% OUTPUTS:
% gw --- Gaussian weights (column vector)

% Last Modified: 07/02/2020 by Xu Zhang
%%
gw = zeros(ng,1);

if ng == 1 % accurate up to polynomial degree p = 1
    
    gw = 1;
    
elseif ng == 4 % accurate up to polynomial degree p = 2
    
    gw = 1/4*ones(4,1);
    
elseif ng == 5 % accurate up to polynomial degree p = 3 
    
    gw(1) = -0.8;
    gw(2:5) = 0.45*ones(4,1);
    
elseif ng == 11 % accurate up to polynomial degree p = 4 
    
    gw(1) = -0.013155555555556;
    gw(2:5) = 0.007622222222222*ones(4,1);
    gw(6:11)= 0.024888888888889*ones(6,1);

elseif ng == 15 % accurate up to polynomial degree p = 5
    
    gw(1) = 0.030283678097089;
    gw(2:5) = 0.006026785714286*ones(4,1);
    gw(6:9) = 0.011645249086029*ones(4,1);
    gw(10:15) = 0.010949141561386*ones(6,1);
    
end