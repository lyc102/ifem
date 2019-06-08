function [Pro,FreeNodec] = transferP2red(elemc,elemf,FreeNode)
%% TRANSFERP2RED
%
%
% Author: Jie Zhou,Jan,17,2013. Discussed with Long Chen

%% Data structure.
if ~exist('FreeNode','var'), FreeNode = []; end 
[elem2dofc,edgec] = dofP2(elemc);
[elem2doff,edgef] = dofP2(elemf);

Nc = max(elemc(:));  Nf = max(elemf(:)); 
Ndofc = Nc + size(edgec,1);Ndoff = Nf + size(edgef,1);
NTc = size(elemc,1);  %NTc:coarse elements number.
%NTf = size(elemf,1);

elem2node2f = sparse(NTc,9); %very important, map from coarse space to fine space
%we just consider the middle points in fine edges. 
elem2node2f(:,1) = elem2doff(1:NTc,4);
elem2node2f(:,2) = elem2doff(1:NTc,5);
elem2node2f(:,3) = elem2doff(1:NTc,6);

elem2node2f(:,4) = elem2doff(NTc+1:2*NTc,4);
elem2node2f(:,5) = elem2doff(NTc+1:2*NTc,5);
elem2node2f(:,6) = elem2doff(NTc+1:2*NTc,6);

elem2node2f(:,7) = elem2doff(2*NTc+1:3*NTc,4);
elem2node2f(:,8) = elem2doff(2*NTc+1:3*NTc,5);
elem2node2f(:,9) = elem2doff(2*NTc+1:3*NTc,6);

%% Assembel the prolongation matrix
%         1      2      3       4       5       6      7       8      9
locRij = [ 0    3/8    3/8      0     -1/8    -1/8      0    -1/8   -1/8 ;...
        -1/8     0    -1/8     3/8     0       3/8   -1/8      0    -1/8 ;...
        -1/8   -1/8     0     -1/8    -1/8      0     3/8     3/8     0  ;...
         1/4     0      0      3/4     1/2      0     3/4      0     1/2 ; ...  
         1/2    3/4     0       0      1/4      0      0      3/4    1/2 ; ...
         1/2    0     3/4       0      1/2     3/4     0       0     1/4];
     
     
% Corase grid order          fine grid order
% 3                          * 
% * *                        * 7
% *   *                      8  *
% *     *                    *    *
% 5       4                  ** 9 **
% *         *                * *    * *
% *           *              2  1   5  4
% *             *            *    * *    *
% 1*****6********2           ** 3 **** 6 **

ii = zeros(40*NTc,1); jj = zeros(40*NTc,1); ss = zeros(40*NTc,1);
index = 0;
for i = 1:6  
    for j = 1:9
        if locRij(i,j)~=0
            ii(index+1:index+NTc) =  elem2dofc(:,i);
            jj(index+1:index+NTc) =  elem2node2f(:,j);
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
ii(ii==0)=[];
jj(jj==0)=[];
ss(ss==0)=[];

% Modification for boundary middle points.
s = accumarray([elem2doff(:,4);elem2doff(:,5);elem2doff(:,6)],1,[Ndoff 1]);
bdEdgeidxf = (s == 1);
bdEdgeidxfinjj = bdEdgeidxf(jj);
ss(bdEdgeidxfinjj) = 2*ss(bdEdgeidxfinjj);
% The transfer operator for 1-6 nodes in the coarse grid is an identical matrix
ii = [(1:Ndofc)'; ii];
jj = [(1:Ndofc)'; jj];
ss = [ones(Ndofc,1); ss];
Pro = sparse(jj,ii,ss,Ndoff,Ndofc);

%% Truncated to free nodes
if ~isempty(FreeNode)
    FreeNodec = FreeNode(FreeNode<=Ndofc);
    Pro = Pro(FreeNode,FreeNodec);
end