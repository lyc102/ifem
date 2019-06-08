function T = auxstructure3(elem)
%% AUXSTRUCTURE3 auxiliary structure for a 3-D triangulation.
%
%  T = AUXSTRUCTURE3(elem) constucts the indices map between elements, 
%  faces, edges and nodes, and the boundary edge information. T is a struct 
%  data. 
%  
%  T.neighor(1:NT,1:4): the indices map of neighbor information of elements, 
%  where neighbor(t,i) is the global index of the element oppoiste to the
%  i-th vertex of the t-th element.
%
%  T.elem2face(1:NT,1:4): the indices map from elements to faces, where 
%  elem2face(t,i) is the global index of the  edge opposite to the i-th 
%  vertex of t-th elemment.
%  
%  T.face(1:NF,1:3): faces, where face(k,i) is the global index of the i-th
%  vetex of the k-th face, and face(k,1)<face(k,2)<face(k,3).
%
%  T.bdFace(1:Nbd,1:3): boundary faces with positive oritentation, where
%  bdFace(k,i) is the global index of the i-th vetex of the k-th boundary
%  face. The positive oritentation means that the order of the three
%  vertices of bdFace(k,:) satisfy the right-handed system and normal
%  points the outside of the domain. Note that this requires elem is
%  positive ordered, i.e., the signed volume of each tetrahedron is
%  positive. If not, use elem = fixorder3(node,elem) to fix the order.
%
%  T.face2elem(1:NF,1:4): the indices map from faces to elements, where 
%  face2elem(k,1:2) are two global indices of the elements sharing the k-th
%  face, and face2elem(k,3:4) are local indices.
%
%  To save space all the data type in T is uint32. When use them as input
%  of sparse(i,j,s,m,n), please change them into double type.
% 
%  See also auxstructure.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

NT = size(elem,1); %N = max(elem(:));
totalFace = uint32([elem(:,[2 3 4]); elem(:,[1 4 3]); ...
                    elem(:,[1 2 4]); elem(:,[1 3 2])]);
totalFace = sort(totalFace,2);
[face, i2, j] = myunique(totalFace);
i1(j(4*NT:-1:1)) = 4*NT:-1:1; i1 = i1';
k1 = ceil(i1/NT); t1 = i1 - NT*(k1-1);
k2 = ceil(i2/NT); t2 = i2 - NT*(k2-1);
idx = (i1 ~= i2); 
neighbor = uint32(accumarray([[t1(idx),k1(idx)];[t2,k2]],[t2(idx);t1],[NT 4]));
elem2face = uint32(reshape(j,NT,4));
face2elem = uint32([t1,t2,k1,k2]);
bdElem = t1(t1 == t2);
bdk1 = k1(t1 == t2);
% consistent ordering of faces is used.
bdFace = [elem(bdElem(bdk1==1),[2 3 4]); elem(bdElem(bdk1==2),[1 4 3]);...
          elem(bdElem(bdk1==3),[1 2 4]); elem(bdElem(bdk1==4),[1 3 2])];
bdFace2elem = [bdElem(bdk1==1);bdElem(bdk1==2);bdElem(bdk1==3);bdElem(bdk1==4)];
T = struct('neighbor',neighbor,'elem2face',elem2face,'face',uint32(face),...
           'face2elem',face2elem,'bdElem',bdElem,'bdFace',uint32(bdFace),'bdFace2elem', bdFace2elem);