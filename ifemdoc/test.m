% A mesh with two tetrahedron with the ascend ordering
elem = [1 4 5 8; 1 4 5 7];
node = [1,0,0; 1,1,1; 1,-1,-1; 0,1,0; -2,-1,0; 1,1,-1; 0,1,1; 0,-1,-1];
NT = size(elem,1);
showmesh3(node,elem,[],'FaceAlpha',0.25);
findelem3(node,elem);
findnode3(node,elem(:));
display(elem);
% generate edge array
totalEdge = uint32([elem(:,[1 2]); elem(:,[1 3]); elem(:,[1 4]); ...
                    elem(:,[2 3]); elem(:,[2 4]); elem(:,[3 4])]);
sortedTotalEdge = sort(totalEdge,2);
[edge, ~, je] = unique(sortedTotalEdge,'rows');
display(edge);
% generate face array
totalFace = uint32([elem(:,[2 3 4]); elem(:,[1 4 3]); ...
                    elem(:,[1 2 4]); elem(:,[1 3 2])]);
sortedTotalFace = sort(totalFace,2);                
[face, ~, jf] = unique(sortedTotalFace,'rows');
display(face);
% generate pointers of indices
elem2edge = uint32(reshape(je,NT,6))
elem2face = uint32(reshape(jf,NT,4))
% find orientation of elem
[v,elemSign] = simplexvolume(node,elem)