function [nodeRho, elemRho, cur, area, v] = curvaturedensity(node,elem, normal, gamma)


n = normal(elem(:,1),:) + normal(elem(:,2),:) + normal(elem(:,3),:);

v12 = node(elem(:,2),:) - node(elem(:,1),:);
v13 = node(elem(:,3),:) - node(elem(:,1),:);

v = cross(v12,v13,2);

area = 0.5*sqrt(sum(v.^2,2));

v = v./[2*area,2*area,2*area];

elemRho = abs(9 - sum(n.^2,2))./area;
elemRho = elemRho +sqrt(eps);


for i = 1:12
    nodeRho = accumarray(elem(:),repmat(elemRho.*area,3,1),[size(node,1),1]);
%     showsolution(node,elem,nodeRho);
%     axis equal;
%     colorbar;
%     pause(0.1);
    b = accumarray(elem(:),[area;area;area],[size(node,1),1]);
    nodeRho =nodeRho./b;
    elemRho = (nodeRho(elem(:,1)) + nodeRho(elem(:,2)) + nodeRho(elem(:,3)))/3;
end

cur = nodeRho;
elemRho = (elemRho/max(elemRho)).^(gamma);
nodeRho = (nodeRho/max(nodeRho)).^(gamma);