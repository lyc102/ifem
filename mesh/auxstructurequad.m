function T = auxstructurequad(elem)
%% AUXSTRUCTUREQUAD auxiliary structure for a 2-D quad mesh
%
%  T = AUXSTRUCTURE(elem) constucts the indices map between elements, edges 
%  and nodes, and the boundary information. T is a structure. 
%
%  T.neighbor(1:NT,1:4): the indices map of neighbor information of elements, 
%  where neighbor(t,i) is the global index of the element oppoiste to the 
%  i-th vertex of the t-th element. 
%
%  T.elem2edge(1:NT,1:4): the indices map from elements to edges, elem2edge(t,i) 
%  is the edge opposite to the i-th vertex of the t-th element.
%
%  T.edge(1:NE,1:2): all edges, where edge(e,i) is the global index of the 
%  i-th vertex of the e-th edge, and edge(e,1) < edge(e,2) 
%
%  T.bdEdge(1:Nbd,1:2): boundary edges with positive oritentation, where
%  bdEdge(e,i) is the global index of the i-th vertex of the e-th edge for
%  i=1,2. The positive oritentation means that the interior of the domain
%  is on the left moving from bdEdge(e,1) to bdEdge(e,2). Note that this
%  requires elem is positive ordered, i.e., the signed area of each
%  triangle is positive. If not, use elem = fixorder(node,elem) to fix the
%  order.
%
%  T.edge2elem(1:NE,1:4): the indices map from edge to element, where 
%  edge2elem(e,1:2) are the global indexes of two elements sharing the e-th
%  edge, and edge2elem(e,3:4) are the local indices of e to edge2elem(e,1:2).
%
%  To save space all the data type in T is uint32. When use them as a input
%  of sparse(i,j,s,m,n), please change them into double type.
%
%  See auxstructure, auxstructure3.
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[NT, NV] = size(elem);
totalEdge = zeros(NV*NT,2);
totalEdge(:,1) = elem(:);
totalEdge(:,2) = elem([NT+1:NV*NT,1:NT]');
totalEdge(:) = sort(totalEdge,2);
matlabversion = version;
if str2double(matlabversion(end-5:end-2)) > 2012
    [edge, i2, j] = unique(totalEdge,'rows','legacy');
else
    [edge, i2, j] = unique(totalEdge,'rows');
end
elem2edge = uint32(reshape(j,NT,NV));
i1(j(NV*NT:-1:1)) = NV*NT:-1:1; 
i1 = i1';
k1 = ceil(i1/NT); 
k2 = ceil(i2/NT); 
t1 = i1 - NT*(k1-1);
t2 = i2 - NT*(k2-1);
ix = (i1 ~= i2); 
neighbor = uint32(accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT NV]));
edge2elem = uint32([t1,t2,k1,k2]);
bdElem = t1(t1 == t2);
bdk1 = k1(t1 == t2);
bdEdge = [elem(bdElem(bdk1==1),[1 2]); elem(bdElem(bdk1==2),[2 3]);...
          elem(bdElem(bdk1==3),[3 4]); elem(bdElem(bdk1==4),[4 1])];
T = struct('neighbor',neighbor,'elem2edge',elem2edge,'edge',edge,...
           'edge2elem',edge2elem,'bdElem',bdElem,'bdEdge',bdEdge);