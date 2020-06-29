function c = mycross(a,b,dim)
%% MYCROSS cross product
%
% c = mycross(a,b) computes the cross product of vectors a and b. It is
% used to replace the built-in cross function which is two times slower.
%
% c = mycross(a,b,dim) the input a, b could be matrix of vectors. The input
% dim specifies the dimension. By default, dim = 2. Namely each row of a
% and b is a vector in R3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('dim','var'), dim = 2; end

% Check dimensions
if (size(a,dim)~=3) || (size(b,dim)~=3)
  error(message('MATLAB:cross:InvalidDimAorBForCrossProd'))
end

c = zeros(size(a,1),3);
c(:,1) = a(:,2).*b(:,3) - a(:,3).*b(:,2);
c(:,2) = a(:,3).*b(:,1) - a(:,1).*b(:,3);
c(:,3) = a(:,1).*b(:,2) - a(:,2).*b(:,1);