function intersectPoint = findintersectbisect(phi,a,b)
% INTERSECTPOINT find the intersect point edge [a,b] and interface phi.
%
% Here use bisection mehtod to find the intersect point.
%
% Huayi Wei <huayiwei1984@gmail.com> 

intersectPoint = (b + a)/2;
phia = phi(a);
phib = phi(b);
phiI = phi(intersectPoint);
isA = false(size(a,1),1);
isB = false(size(a,1),1);

v = b - a; 
h = sqrt(sum(v.*v,2));
epsh = eps.*h;
isNotOK = h > epsh & phiI ~= 0;

while any(isNotOK)

intersectPoint(isNotOK,:) = (b(isNotOK,:)+a(isNotOK,:))/2;
phiI(isNotOK) = phi(intersectPoint(isNotOK,:));


isA(isNotOK) = phia(isNotOK).*phiI(isNotOK) > 0;
isB(isNotOK) = phib(isNotOK).*phiI(isNotOK) > 0;

a(isA,:) = intersectPoint(isA,:);
b(isB,:) = intersectPoint(isB,:);

phia(isA) = phiI(isA);
phib(isB) = phiI(isB);

h(isNotOK) = h(isNotOK)/2;

isNotOK(isNotOK) = h(isNotOK) > epsh(isNotOK) & phiI(isNotOK) ~= 0;

isA(:) = false;
isB(:) = false;
end