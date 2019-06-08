function [normal,edgeLength,unitNormal] = edgenormal(node,edge)
% compute normal vector of a triangular face

edgeVec = node(edge(:,2),:)-node(edge(:,1),:);
normal = [edgeVec(:,2), -edgeVec(:,1)];
if nargout > 1
    edgeLength = sqrt(sum(normal.^2,2));
    if nargout > 2
        unitNormal = normal./repmat(edgeLength,[1,2]);
    end
end