function plotcircle(x0, y0, r)
theta = 0:2*pi/1000:2*pi;
x = x0 + r*cos(theta);
y = y0 + r*sin(theta);
plot(x, y,'r-','LineWidth',2);