function q = quadelem2node(node,elem,func1,func2)

%q(:,i) = \int_K f1 f2 \lamba_i dx
%where i is the local index of the node

NT = size(elem,1);
b_x=[1/3;(6+sqrt(15))/21;(9-2*sqrt(15))/21;(6+sqrt(15))/21
     (6-sqrt(15))/21;(9+2*sqrt(15))/21;(6-sqrt(15))/21];
b_y=[1/3;(6+sqrt(15))/21;(6+sqrt(15))/21;(9-2*sqrt(15))/21
        (6-sqrt(15))/21;(6-sqrt(15))/21;(9+2*sqrt(15))/21];
w=[9/80;(155+15^(1/2))/2400;(155-15^(1/2))/2400];

n1 = elem(:,1);	 n2 = elem(:,2);   n3 = elem(:,3);
x1 = node(n1,1); y1 = node(n1,2);   
x2 = node(n2,1); y2 = node(n2,2);    
x3 = node(n3,1); y3 = node(n3,2);

% aa(:,1) = y2-y3;	 bb(:,1) = x3-x2;	
% aa(:,2) = y3-y1;	 bb(:,2) = x1-x3;	
% aa(:,3) = y1-y2;	 bb(:,3) = x2-x1;

%the inverse of B where \hat{x} = (B^{-1})^T(x - c)
%the 3rd dimension follows the index of elem
%Binv(1,:,:) = aa(:,1:2)';
%Binv(2,:,:) = bb(:,1:2)';
%Binv(1,1,:) = aa(:,1); Binv(1,2,:) = aa(:,2);
%Binv(2,1,:) = bb(:,1); Binv(2,2,:) = bb(:,2);

area = (-x2.*y1+x3.*y1+x1.*y2-x3.*y2-x1.*y3+x2.*y3)/2;

xx = x1*b_x'+x2*b_y'+x3*(1-b_x-b_y)';
yy = y1*b_x'+y2*b_y'+y3*(1-b_x-b_y)';

%shape function on the reference element
%lambda1 = @(p)p(:,1);
%lambda2 = @(p)p(:,2);
%lambda3 = @(p)1 - p(:,1) - p(:,2);

if nargin==4
    func = @(p) func1(p).*func2(p);
else
    func = @(p) func1(p);
end

f = func([xx(:) yy(:)]);
f = reshape(f,size(xx));
%f = eval(func);
f1 = f.*repmat(b_x',NT,1);
f2 = f.*repmat(b_y',NT,1);
f3 = f.*repmat((1-b_x-b_y)',NT,1);

q(:,1) = 2*area.*(w(1)*f1(:,1)+w(2)*(f1(:,2)+f1(:,3)+f1(:,4))...
    +w(3)*(f1(:,5)+f1(:,6)+f1(:,7))); 
q(:,2) = 2*area.*(w(1)*f2(:,1)+w(2)*(f2(:,2)+f2(:,3)+f2(:,4))...
    +w(3)*(f2(:,5)+f2(:,6)+f2(:,7)));
q(:,3) = 2*area.*(w(1)*f3(:,1)+w(2)*(f3(:,2)+f3(:,3)+f3(:,4))...
    +w(3)*(f3(:,5)+f3(:,6)+f3(:,7)));