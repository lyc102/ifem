function z = phi1(p,a) % level set function
anu = 0.02*sqrt(5);
theta = atan2((p(:,2)-anu),(p(:,1)-anu));
z = (p(:,2)-anu).^2 + (p(:,1)-anu).^2 - (0.5+0.2*sin(a*theta)).^2;
end
