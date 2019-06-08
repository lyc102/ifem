function [Pro,freeEdgec] = transferCRred(elemc,elemf,freeEdge)
%% TRANSFERCRRED prolongation operator of CR element
% 
% Pro = TRANSFERCRRED(elemc,elemf) constructs a prolongation operator of CR
% linear elements from the coarse grid elemc to the fine grid elemf. The
% fine grid is a uniform (red) refinement of the coarse grid. Note that
% since the CR element spaces are not nested, the prolongation operator for
% C-R element is not unique.
%
% [Pro,freeEdgec] = TRANSFERCRRED(elemc,elemf,freeEdge) restricts the
% prolongation operator to free edges. 
%
% Author: Jie Zhou and Long Chen. m-lint: Long Chen.

% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if ~exist('freeEdge','var'), freeEdge = []; end 

%% Data structure
% edge and elem2edge
[elem2edgec,edgec] = dofedge(elemc);
[elem2edgef,edgef] = dofedge(elemf);
NEf = size(edgef,1);
NEc = size(edgec,1);   
NTc = size(elemc,1);
Nc  = max(elemc(:));
% index map inside one element
elem2dofc2f = zeros(NTc,9);
elem2dofc2f(:,1) = elem2edgef(1:NTc,1);
elem2dofc2f(:,2) = elem2edgef(1:NTc,2);
elem2dofc2f(:,3) = elem2edgef(1:NTc,3);
elem2dofc2f(:,4) = elem2edgef(NTc+1:2*NTc,1);
elem2dofc2f(:,5) = elem2edgef(NTc+1:2*NTc,2);
elem2dofc2f(:,6) = elem2edgef(NTc+1:2*NTc,3);
elem2dofc2f(:,7) = elem2edgef(2*NTc+1:3*NTc,1);
elem2dofc2f(:,8) = elem2edgef(2*NTc+1:3*NTc,2);
elem2dofc2f(:,9) = elem2edgef(2*NTc+1:3*NTc,3);

% corase grid index         fine grid index
% *                          * 
% * *                        * 7
% *   *                      8  *
% *     *                    * T3  *
% 2       1                  ** 9  **
% *         *                * *    * *
% *           *              2  1   5   4
% *             *            * T1 * * T2  *
% ***** 3 ********           ** 3 **** 6  **    

%% Local restriction operator
%          1     2     3     4     5      6     7     8     9           
locRij = [ 0   -0.5  -0.5    1    0.5    0.5    1    0.5   0.5 ;...
          0.5    1    0.5  -0.5    0    -0.5   0.5    1    0.5 ;...
          0.5   0.5    1    0.5   0.5     1   -0.5  -0.5    0 ];    
 
%% Assemble the prolongation matrix
ii = zeros(24*NTc,1); jj = zeros(24*NTc,1); ss = zeros(24*NTc,1);
index = 0;
for i = 1:3  
    for j = 1:9
        if locRij(i,j)~=0
            ii(index+1:index+NTc) = elem2edgec(:,i);
            jj(index+1:index+NTc) = elem2dofc2f(:,j);
            if (j == 1) || (j == 5) || (j == 9)
                ss(index+1:index+NTc) = locRij(i,j);
            else
                ss(index+1:index+NTc) = locRij(i,j)/2;
            end
%         Every interior edges of the coarse grid will be computed twice.
%         So the corresponding weight is halved. Values for boundary edges
%         will be modified later.
             index = index + NTc;
        end
    end
end

% Double the value for boundary edges
s = accumarray(elem2edgef(:), 1, [NEf 1]);
bdEdgeidxf = (s == 1);
bdEdgeidxfinjj = bdEdgeidxf(jj);
ss(bdEdgeidxfinjj) = 2*ss(bdEdgeidxfinjj);
Pro = sparse(double(jj),double(ii),ss,NEf,NEc);

%% Truncated to free edges
if ~isempty(freeEdge)
    isFixedEdge  = true(NEf,1);
    isFixedEdge(freeEdge) = false;
    fixedEdgeNode  = edgef(isFixedEdge,:);
    isFixedEdgec = false(NEc,1);
    index = fixedEdgeNode(:)-Nc;
    index(index<=0) = [];
    isFixedEdgec(index) = true;
    freeEdgec = find(~isFixedEdgec);
    Pro = Pro(freeEdge,freeEdgec);
end