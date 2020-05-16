function [flag, ixErrElem] = checkpoly(node,elem)
%% CHECKPOLY checks if a polygonal mesh is correctly constructed
% flag = # of polygons with wrong indexing
% TODO: vectorize using nested cellfun
%
if ~iscell(elem); elem = num2cell(elem,2); end
NT = size(elem,1);
flag = false(NT,1);

for i = 1:NT
    [~, isOnElem] = inpolygon(node(:,1),node(:,2),node(elem{i},1),node(elem{i},2));
    flag(i) =  ~isempty(setdiff(find(isOnElem), elem{i}));
end
flag = sum(flag);
ixErrElem = find(flag);

end