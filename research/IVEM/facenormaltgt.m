function [normal,area,unitNormal,unittgt1,unittgt2] = facenormaltgt(node,face)
% compute normal vector of a triangular face

v12 = node(face(:,2),:)-node(face(:,1),:);
v13 = node(face(:,3),:)-node(face(:,1),:);
normal = mycross(v12,v13,2);
if nargout > 1
    area = sqrt(sum(normal.^2,2))/2;
    if nargout > 2
        unitNormal = 0.5*normal./repmat(area,[1,3]);
        if nargout > 3
            unittgt1 = v12./(sum(v12.^2,2).^(1/2)*ones(1,3));
            unittgt2 = mycross(unitNormal,unittgt1,2);
        end
    end
end