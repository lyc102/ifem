function n = intlambda(alpha,d,volume)
%% INTLAMBDA integral of barycentric coordinate on a simplex
%
% n = INTLAMBDA(alpha,d) computes the following integral 
%
%   int_T lambda^alpha dx = alpha!d!/(|alpha|+d)! |T|
% - T: a simplex, assuming |T| = 1 if the third argument is None.
% - alpha: a multi-index counter that allows repetition
%          e.g. [1 1 3] results a multi-index [2 0 1]'
% - d: dimension of the underlying space
%      d=1, line; d=2, triangle; d=3, tetrahedron.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if isrow(alpha); alpha = alpha'; end % convert to column vector if row

alpha = accumarray(alpha,1); % convert alpha to a multi-index
n = prod(factorial(alpha))*factorial(d)/factorial(sum(alpha)+d);
if nargin >= 3; n = n*volume; end