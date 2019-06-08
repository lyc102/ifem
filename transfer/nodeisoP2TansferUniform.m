function [Pro,FreeNodec] = nodeisoP2TansferUniform(elemc,elemf,FreeNode)
%  Only works for uniformrefinered.
%
%
% Created by Jie Zhou,Jan,17,2013.  
% Talk with Long Chen
%% Data structure.
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
elem2node2f(:,1)=elem2doff(1:NTc,4);
elem2node2f(:,2)=elem2doff(1:NTc,5);
elem2node2f(:,3)=elem2doff(1:NTc,6);

elem2node2f(:,4)=elem2doff(NTc+1:2*NTc,4);
elem2node2f(:,5)=elem2doff(NTc+1:2*NTc,5);
elem2node2f(:,6)=elem2doff(NTc+1:2*NTc,6);

elem2node2f(:,7)=elem2doff(2*NTc+1:3*NTc,4);
elem2node2f(:,8)=elem2doff(2*NTc+1:3*NTc,5);
elem2node2f(:,9)=elem2doff(2*NTc+1:3*NTc,6);

%% Assembel the transfer matrix
%        1       2         3        4       5       6         7           8          9
c0=1/2;
locPij=[ 0    1/2*c0    1/2*c0      0        0        0          0         0        0;...
         0    0         0         1/2*c0     0       1/2*c0      0         0        0;...
         0    0         0           0        0        0         1/2*c0     1/2*c0   0;...
         0    0         0         1/2*c0     1/2      0         1/2*c0     0        1/2; ...  
         1/2  1/2*c0    0           0        0        0          0         1/2*c0   1/2; ...
         1/2  0         1/2*c0      0        1/2      1/2*c0     0         0        0];
%%Corase grid order          fine grid order
% 3                          * 
% * *                        *  7
% *   *                      8    *
% *     *                    *      *
% 5       4                  ***9******
% *         *                *      *    4
% *           *              2  1   5      *
% *             *            *      *        *
% 1*****6********2           ***3********6*****     

ii = zeros(40*NTc,1); jj = zeros(40*NTc,1); ss = zeros(40*NTc,1);
index = 0;
for i = 1:6  
    for j = 1:9
        if locPij(i,j)~=0
            %Pij = locPij(i,j)*elem2edgeSignc(:,i).*elem2edgesignc2f(:,j);
            ii(index+1:index+NTc) =  elem2dofc(:,i);
            jj(index+1:index+NTc) =  elem2node2f(:,j);
            ss(index+1:index+NTc) =  locPij(i,j);
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
ii = [(1:Ndofc)'; ii];
jj = [(1:Ndofc)'; jj];
ss = [ones(Ndofc,1); ss];
Pro     =  sparse(jj,ii,ss,Ndoff,Ndofc);
% Pro(1:Ndofc,1:Ndofc) = speye(Ndofc,Ndofc);
% The transfer for  is a identical matrix, so the local matrix is 6 by 9.


%% Truncated to free nodes
if ~isempty(FreeNode)
    FreeNodec = FreeNode(FreeNode<=Ndofc);
    Pro = Pro(FreeNode,FreeNodec);
end