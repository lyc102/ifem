function T = auxstructurehex(elem)
%% AUXSTRUCTUREHEX auxiliary structure for a 3-D HEX mesh.
%
%  T = AUXSTRUCTUREHEX(elem) constucts the indices map between elements, 
%  faces, edges and nodes, and the boundary edge information. T is a struct 
%  data. 
%  
%  T.neighor(1:NT,1:6): the indices map of neighbor information of elements, 
%  where neighbor(t,i) is the global index of the element oppoiste to the
%  i-th vertex of the t-th element.
%
%  T.elem2face(1:NT,1:6): the indices map from elements to faces, where 
%  elem2face(t,i) is the global index of the  edge opposite to the i-th 
%  vertex of t-th elemment.
%  
%  T.face(1:NF,1:4): faces, where face(k,i) is the global index of the i-th
%  vetex of the k-th face, and face(k,1)<face(k,2)<face(k,3)<face(k,4).
%
%  T.bdFace(1:Nbd,1:4): boundary faces with positive oritentation, where
%  bdFace(k,i) is the global index of the i-th vetex of the k-th boundary
%  face. The positive oritentation means that the order of the three
%  vertices of bdFace(k,:) satisfy the right-handed system and normal
%  points the outside of the domain. 
%
%  T.face2elem(1:NF,1:4): the indices map from faces to elements, where 
%  face2elem(k,1:2) are two global indices of the elements sharing the k-th
%  face, and face2elem(k,3:4) are local indices of e to edge2elem(k,1:2).
%
%  To save space all the data type in T is uint32. When use them as a input
%  of sparse(i,j,s,m,n), please change them into double type.
% 
%  See also auxstructure, auxstructurequad, auxstructurehex.
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

NT = size(elem,1);
totalFace = uint32([elem(:,[1 4 3 2]);elem(:,[1 2 6 5]);elem(:,[5 6 7 8]);elem(:,[8 7 3 4]);...
               elem(:,[4 1 5 8]); elem(:,[2 3 7 6])]);
matlabversion = version;
if str2double(matlabversion(end-5:end-2)) > 2012
    [face, i2, j] = unique(sort(totalFace,2),'rows','legacy');
else           
    [face, i2, j] = unique(sort(totalFace,2),'rows');
end
i1(j(6*NT:-1:1)) = 6*NT:-1:1; i1 = i1';
k1 = ceil(i1/NT); t1 = i1 - NT*(k1-1);
k2 = ceil(i2/NT); t2 = i2 - NT*(k2-1);
idx = (i1 ~= i2); 
neighbor = uint32(accumarray([[t1(idx),k1(idx)];[t2,k2]],[t2(idx);t1],[NT 6]));
elem2face = uint32(reshape(j,NT,6));
face2elem = uint32([t1,t2,k1,k2]);
bdElem = t1(t1 == t2);
bdk1 = k1(t1 == t2);
% consistent ordering of faces is used.
bdFace = [elem(bdElem(bdk1==1),[1 4 3 2]); elem(bdElem(bdk1==2),[1 2 6 5]);...
          elem(bdElem(bdk1==3),[5 6 7 8]); elem(bdElem(bdk1==4),[8 7 3 4]);...
          elem(bdElem(bdk1==5),[4 1 5 8]); elem(bdElem(bdk1==6),[2 3 7 6])];
bdFace2elem = [bdElem(bdk1==1)
               bdElem(bdk1==2)
               bdElem(bdk1==3)
               bdElem(bdk1==4)
               bdElem(bdk1==5)
               bdElem(bdk1==6)];
T = struct('neighbor',neighbor,'elem2face',elem2face,'face',uint32(face),...
           'face2elem',face2elem,'bdElem',bdElem,'bdFace',uint32(bdFace),...
           'bdFace2elem', bdFace2elem);