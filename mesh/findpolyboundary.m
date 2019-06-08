function [bdNode,bdFace,isBdNode] = findpolyboundary(face)

NF = size(face, 1);
[uface, i2, j] = myunique(sort(face, 2)); %#ok<*ASGLU>
i1(j(NF:-1:1)) = NF:-1:1;
i1 = i1';
bdFace = face(i1(i1 == i2), :);

isTriFace = bdFace(:, 4) == 0;
isBdNode(bdFace(isTriFace,1:3)) = true;
isBdNode(bdFace(~isTriFace, :)) = true;
bdNode = find(isBdNode);
