function z = ellipse(p)

x = p(:,1);
y = p(:,2);

z = x.^2+3*y.^2 - 0.5;