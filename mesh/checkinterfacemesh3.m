function flag = checkinterfacemesh3(node, elem, interfaceData)
%% CHECKINTERFACEMESH3 check the validity of the polyhedra mesh generated
%   by interfacemesh3.
%
%   Here we check  validity of the polyhedra mesh by the Euler formla: 
%            F - E + V = 2 
%   where F, E and V are the number of faces, edges and vertices of a simply
%   closed polyhedron.
%
%   Author: Huayi Wei <weihuayi@xtu.edu.cn>, based on discussion with Long Chen.
%
%   See also: interfacemesh3
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

face = interfaceData.face;
face2elem = interfaceData.face2elem;
N = size(node,1);
NF = size(face,1);

NP = max(face2elem);
isCutPoly = false(NP, 1);
isCutPoly(face2elem) = true;
cutPolyIdx = zeros(NP, 1);
NP = sum(isCutPoly);
cutPolyIdx(isCutPoly) = 1:NP;
face2elem = cutPolyIdx(face2elem);

isTriFace = face(:,4) == 0;
edge = 4*ones(NF,1);
edge(isTriFace) = 3;
E = accumarray(face2elem, edge)/2;
F = accumarray(face2elem, ones(NF,1));

poly2node = sparse(face2elem(isTriFace)*ones(1,3), face(isTriFace, 1:3), 1, NP,N)...
+ sparse(face2elem(~isTriFace)*ones(1,4), face(~isTriFace, 1:4), 1, NP,N);
poly2node = poly2node > 0;
NV = poly2node*ones(N,1);

isBdPoly = ((F - E + NV) ~= 2);

if ~isempty(find(isBdPoly,1))
    figure
    showmesh(node, face(isTriFace & isBdPoly(face2elem),1:3));
    sface = face(~isTriFace & isBdPoly(face2elem),1:4);
    hold on
    tsface = [sface(:, 1:3);sface(:,[3, 4,1])];
    showmesh(node, tsface, 'Facecolor','y');
    
%     bdPolyIdx = find(isBdPoly);
%     figure
%     for i = 1:size(bdPolyIdx)
%         sface = face(~isTriFace & face2elem == bdPolyIdx(i),1:4);
%         tsface = [sface(:, 1:3);sface(:,[3, 4,1])];
%         showmesh(node, [face(isTriFace & face2elem == bdPolyIdx(i),1:3);tsface]);
%     end
    flag = false;
else
    flag = true;
end




