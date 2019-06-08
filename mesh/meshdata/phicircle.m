function z = circle(p)

x = p(:,1);
y = p(:,2);
z = sqrt((x- 0.5).^2+y.^2) - 0.01;
