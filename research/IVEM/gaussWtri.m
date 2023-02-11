function gw = gaussWtri(ng)

%% USAGE: generate Gaussian weights on a triangle
%
% INPUTS:
% nguass --- number of Gaussian nodes in an triangle
%              possible values = 1, 3, 4, 6, 7, 12, 13.
%
% OUTPUTS:
% gw --- Gaussian weights (column vector)

% Last Modified: 07/05/2020 by Xu Zhang

%%
gw = zeros(ng,1);

if ng == 1 % accurate up to polynomial degree p = 1
    
    gw = 1;
    
elseif ng == 3 % accurate up to polynomial degree p = 2
    
    gw(1:3) = 1/3; 
    
elseif ng == 4 % accurate up to polynomial degree p = 3 
    
    gw(1)  = -0.5625;
    gw(2:4) = 0.520833333333333;
    
elseif ng == 6 % accurate up to polynomial degree p = 4 
    
    gw(1:3) = 0.109951743655322;
    gw(4:6) = 0.223381589678011;

elseif ng == 7 % accurate up to polynomial degree p = 5
    
    gw(1)   = 9/40;
    gw(2:4) = (155-sqrt(15))/1200;
    gw(5:7) = (155+sqrt(15))/1200;
    
elseif ng == 12 % accurate up to polynomial degree p = 6 
    
    gw(1:3) = 0.050844906370207;
    gw(4:6) = 0.116786275726379;
    gw(7:12)= 0.082851075618374;
    
elseif ng == 13 % accurate up to polynomial degree p = 7 
    
    gw(1)  = -0.149570044467670;
    gw(2:4) = 0.175615257433204;
    gw(5:7) = 0.053347235608839;  
    gw(8:13)= 0.077113760890257;
    
end
