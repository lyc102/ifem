function q = triangleratio(node,elem)
a=sqrt(sum((node(elem(:,2),:)-node(elem(:,1),:)).^2,2));
b=sqrt(sum((node(elem(:,3),:)-node(elem(:,1),:)).^2,2));
c=sqrt(sum((node(elem(:,3),:)-node(elem(:,2),:)).^2,2));
% radius of inscribed circle
% formula: r = sqrt [ (s - a)(s - b)(s - c) / s ] where s = (a + b + c) / 2
r=1/2*sqrt((b+c-a).*(c+a-b).*(a+b-c)./(a+b+c));
% radius of the circumcircle of a triangle
% formula: R = a*b*c/(4*Area) where Area = sqrt [ s (s - a)(s - b)(s - c) ]
R=a.*b.*c./sqrt((a+b+c).*(b+c-a).*(c+a-b).*(a+b-c));
q=2*r./R;
end
