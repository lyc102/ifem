function z = phi5(p)
x = p(:,1); y = p(:,2);
r = x.^2 + y.^2;
theta = atan2(y,x);
isNeg = theta < 0;
theta(isNeg) = theta(isNeg) + 2*pi;

x1 = 16*sin(theta).^3;
y1 = 13*cos(theta) - 5*cos(2*theta) - 2*cos(3*theta) - cos(4*theta);
z = r - (x1.^2 + y1.^2);
end