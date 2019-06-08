function [normal,area,unitNormal] = facenormal(node,face)
% compute normal vector of a triangular face

v12 = node(face(:,2),:)-node(face(:,1),:);
v13 = node(face(:,3),:)-node(face(:,1),:);
normal = mycross(v12,v13,2);
if nargout > 1
    area = sqrt(sum(normal.^2,2))/2;
    if nargout > 2
        unitNormal = 0.5*normal./repmat(area,[1,3]);
    end
end