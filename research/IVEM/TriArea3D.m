function A = TriArea3D(X1,X2,X3)

% Xi n*3 are the three points of triangles where n is the number of
% triangle

x1 = X1(:,1);  y1 = X1(:,2);   z1 = X1(:,3);
x2 = X2(:,1);  y2 = X2(:,2);   z2 = X2(:,3);
x3 = X3(:,1);  y3 = X3(:,2);   z3 = X3(:,3);
A = 1/2*(((x1-x3).*(y2-y1) - (x1-x2).*(y3-y1)).^2 + ...
    ((y1-y3).*(z2-z1) - (y1-y2).*(z3-z1)).^2 + ...
    ((z1-z3).*(x2-x1) - (z1-z2).*(x3-x1)).^2).^(1/2);

A = abs(A);

end