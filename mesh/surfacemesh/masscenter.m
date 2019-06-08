function   [massCenter] = masscenter(node,elem,surfacedata,gamma)
%
N= size(node,1);
b = zeros(N,1);
ax = zeros(N,1);
ay = zeros(N,1);
az = zeros(N,1);

c = tricircumcenter3d(node,elem);
c = updatecircumcenter(node,elem,c);

normal = surfacedata.unitoutnormal(node);
[~, rho] = curvaturedensity(node,elem, normal,gamma);

mc1 = (node(elem(:,2),:) + node(elem(:,3),:))/2;

mc2 = (node(elem(:,3),:) + node(elem(:,1),:))/2;

mc3 = (node(elem(:,1),:) + node(elem(:,2),:))/2;

%%%%% first vertex
tmpbc1 = (mc2 + c + node(elem(:,1),:))/3;
tmpv1 = cross(c - node(elem(:,1),:), mc2 - node(elem(:,1),:),2);
tmpArea1 = 0.5*sqrt(sum(tmpv1.^2,2)).*rho;

tmpbc2 = (mc3 + c + node(elem(:,1),:))/3;
tmpv2 = cross(mc3 - node(elem(:,1),:), c - node(elem(:,1),:), 2);
tmpArea2 = 0.5*sqrt(sum(tmpv2.^2,2)).*rho;


b = b + accumarray(elem(:,1),tmpArea1 + tmpArea2,[N,1]);
ax = ax + accumarray(elem(:,1),tmpbc1(:,1).*tmpArea1 + tmpbc2(:,1).*tmpArea2,[N,1]);
ay = ay + accumarray(elem(:,1),tmpbc1(:,2).*tmpArea1 + tmpbc2(:,2).*tmpArea2,[N,1]);
az = az + accumarray(elem(:,1),tmpbc1(:,3).*tmpArea1 + tmpbc2(:,3).*tmpArea2,[N,1]);

%%%% second vertex
tmpbc1 = (mc3 + c + node(elem(:,2),:))/3;
tmpv1 = cross(c - node(elem(:,2),:), mc3 - node(elem(:,2),:),2);
tmpArea1 = 0.5*sqrt(sum(tmpv1.^2,2)).*rho;

tmpbc2 = (mc1 + c + node(elem(:,2),:))/3;
tmpv2 = cross(mc1 - node(elem(:,2),:), c - node(elem(:,2),:), 2);
tmpArea2 = 0.5*sqrt(sum(tmpv2.^2,2)).*rho;

b = b + accumarray(elem(:,2),tmpArea1 + tmpArea2,[N,1]);
ax = ax + accumarray(elem(:,2),tmpbc1(:,1).*tmpArea1 + tmpbc2(:,1).*tmpArea2,[N,1]);
ay = ay + accumarray(elem(:,2),tmpbc1(:,2).*tmpArea1 + tmpbc2(:,2).*tmpArea2,[N,1]);
az = az + accumarray(elem(:,2),tmpbc1(:,3).*tmpArea1 + tmpbc2(:,3).*tmpArea2,[N,1]);

%%%% third vertex

tmpbc1 = (mc1 + c + node(elem(:,3),:))/3;
tmpv1 = cross(c - node(elem(:,3),:), mc1 - node(elem(:,3),:),2);
tmpArea1 = 0.5*sqrt(sum(tmpv1.^2,2)).*rho;

tmpbc2 = (mc2 + c + node(elem(:,3),:))/3;
tmpv2 = cross(mc2 - node(elem(:,3),:), c - node(elem(:,3),:), 2);
tmpArea2 = 0.5*sqrt(sum(tmpv2.^2,2)).*rho;

b = b + accumarray(elem(:,3),tmpArea1 + tmpArea2,[N,1]);
ax = ax + accumarray(elem(:,3),tmpbc1(:,1).*tmpArea1 + tmpbc2(:,1).*tmpArea2,[N,1]);
ay = ay + accumarray(elem(:,3),tmpbc1(:,2).*tmpArea1 + tmpbc2(:,2).*tmpArea2,[N,1]);
az = az + accumarray(elem(:,3),tmpbc1(:,3).*tmpArea1 + tmpbc2(:,3).*tmpArea2,[N,1]);


massCenter = [ax./b, ay./b,az./b];

l = sum((massCenter - node).*normal,2);
massCenter = massCenter - [l,l,l].*normal;