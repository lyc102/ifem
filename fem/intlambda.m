function n = intlambda(alpha,d)
%% INTLAMBDA integral of barycentric coordinate
%
% n = INTLAMBDA(alpha,d) computes the following integral 
%
%   int_T lambda^alpha dx = alpha!d!/(|alpha|+d)!
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if size(alpha,2) > size(alpha,1) % row vector
    alpha = alpha';              % transfer to column vector
end
alpha = accumarray(alpha,1);
n = prod(factorial(alpha))*factorial(d)/factorial(sum(alpha)+d);