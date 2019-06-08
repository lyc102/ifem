function z = msign(x)
%% MSIGN modified sign function
%
% z = MSIGN(x) returns zero when abs(x)<sqrt(eps)

z = sign(x);
idx = (abs(x)< sqrt(eps));
z(idx) = 0;