function c = updatecircumcenter(node,elem,c)

v1 = node(elem(:,1),:);
v2 = node(elem(:,2),:);
v3 = node(elem(:,3),:);

v21 = v2 - v1;
v31 = v3 - v1;

v32 = v3 - v2;

L21 = sum(v21.^2,2);
L31 = sum(v31.^2,2);
L32 = sum(v32.^2,2);

flag1 = L21 + L31 < L32;
flag2 = L32 + L21 < L31;
flag3 = L31 + L32 < L21;

c(flag1,:) = (v2(flag1,:) + v3(flag1,:))/2;
c(flag2,:) = (v1(flag2,:) + v3(flag2,:))/2;
c(flag3,:) = (v1(flag3,:) + v2(flag3,:))/2;
