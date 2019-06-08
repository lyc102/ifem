function Pro = transferedgebisect(elemc,elemf,tree,bdFlagc,bdFlagf)
%% TRANSFEREDGEBISECT transfer operator between bisected grids
%
%  Works for biscetion grids.
%
%  input elemc,elemf,tree
%  ouput Pro,Res
%   
% Added by Jie Zhou. Based on discussion with Long Chen.
%  
%
%% Data structure
% elem2edge and elem2edgeSign
[elem2edgec,edgec,elem2edgeSignc] = dofedge(elemc);
[elem2edgef,edgef,elem2edgeSignf] = dofedge(elemf);
elem2edgeSignf = double(elem2edgeSignf);
elem2edgeSignc = double(elem2edgeSignc);
% number of edges and elements
NEc = size(edgec,1);  NEf = size(edgef,1); 
NTc = size(elemc,1);  NTf = size(elemf,1);
% index map of left to right
elemcl2r = (1:NTc)'; % a map from left child to the last right child
elemcl2r1 = (1:NTc)'; % a map from left child to the first right child
elemcl2r(tree(:,1)) = double(tree(:,3));  
elemcl2r1(tree(end:-1:1,1)) = double(tree(end:-1:1,3));
% Five cases: only works for bisect not coarsen
% assumption: new nodes are appended to the coarse one. but coarsen may
% remove nodes with smaller index.
isBisect = false(NTf,1);
isBisect(tree(:,1)) = true;
isNotBisect = ~isBisect(1:NTc);
splite0 = find(isNotBisect);  
splite1 = find((elemc(1:NTc,1) == elemf(1:NTc,2)) & ...
               (elemc(1:NTc,1) == elemf(elemcl2r(1:NTc),3)));
splite2 = find((elemc(1:NTc,1) == elemf(1:NTc,2)) & ...
               (elemc(1:NTc,1) ~= elemf(elemcl2r(1:NTc),3)));
splite3 = find((elemc(1:NTc,1) ~= elemf(1:NTc,2)) & ...
               (elemc(1:NTc,1) == elemf(elemcl2r1(1:NTc),3)));
splite4 = find((elemc(1:NTc,1) ~= elemf(1:NTc,2)) & ...
               (elemc(1:NTc,1) ~= elemf(elemcl2r1(1:NTc),3)) & ...
                   isBisect(1:NTc));
if elemcl2r == elemcl2r1  % only one bisection
    splite1 = tree(:,2);
    splite2 = [];  splite3 = [];   splite4 = [];    
end

%% Assemble the prolongation matrix
% preallocation
ii = zeros(14*NTc,1); jj = zeros(14*NTc,1); ss = zeros(14*NTc,1);
index = 0;
%
% Every interior edges of the coarse grid will be computed twice. So in
% locRij theses coarse edge is half of the original value. Values for
% boundary edges will be modified later.

%**************************************************************************
%    CASE ZERO: the coarse triangle is not splitted
%**************************************************************************
NTc0 = length(splite0);
if NTc0 > 0
    elem2edgec2f(splite0,1) = elem2edgef(splite0,1);
    elem2edgec2f(splite0,2) = elem2edgef(splite0,2);
    elem2edgec2f(splite0,3) = elem2edgef(splite0,3);

    locRij = [1   0    0  ;...
              0   1    0  ;...
              0   0    1 ];
    for i = 1:3  
        Rij  = locRij(i,i);
        ii(index+1:index+NTc0) = elem2edgec(splite0,i);
        jj(index+1:index+NTc0) = elem2edgec2f(splite0,i);
        ss(index+1:index+NTc0) = Rij/2;
        index = index + NTc0;
    end
end

%**************************************************************************
%    CASE ONE: one bisection
%**************************************************************************
%                      V1                  V1       ef1:[V1,V2]
%                     * *                  *|*      ef2:[V2,V4]
%                    *   *                * | *     ef3:[V4,V3]
%                   *     *              *  |  *    ef4:[V3,V1]
%                 ec3      ec2         ef1 ef5 ef4  ef5:[V1,V4] 
%                 *  T=L    *          * L  | R  *  
%                *           *        *     |      *
%              V2 *** ec1 *** V3    V2*ef2* V4*ef3*V3  

NTc1 = length(splite1);  
if NTc1 > 0
    elem2edgec2f(splite1,1) = elem2edgef(splite1,1);
    elem2edgec2f(splite1,2) = elem2edgef(splite1,2);
    elem2edgec2f(splite1,3) = elem2edgef(elemcl2r(splite1),3);
    elem2edgec2f(splite1,4) = elem2edgef(elemcl2r(splite1),1);
    elem2edgec2f(splite1,5) = elem2edgef(elemcl2r(splite1),2);

    elem2edgesignc2f(splite1,1) = elem2edgeSignf(splite1,1);
    elem2edgesignc2f(splite1,2) = elem2edgeSignf(splite1,2);
    elem2edgesignc2f(splite1,3) = elem2edgeSignf(elemcl2r(splite1),3);
    elem2edgesignc2f(splite1,4) = elem2edgeSignf(elemcl2r(splite1),1);
    elem2edgesignc2f(splite1,5) = elem2edgeSignf(elemcl2r(splite1),2);
%             1     2     3    4    5
    locRij = [0   0.5   0.5    0    0  ;...
              0     0     0    1   -1 ;...
              1     0     0    0    1 ];
    for i = 1:3  
        for j = 1:5
            if locRij(i,j) ~= 0
                Rij  = locRij(i,j)*elem2edgeSignc(splite1,i).*elem2edgesignc2f(splite1,j);
                ii(index+1:index+NTc1) = elem2edgec(splite1,i);
                jj(index+1:index+NTc1) = elem2edgec2f(splite1,j);
                ss(index+1:index+NTc1) = Rij/2;
                index = index + NTc1;
            end
        end
    end
end
%**************************************************************************
%    CASE TWO: two bisections. right child is bisect.
%**************************************************************************
%                      V1                  V1          ef1:[V2,V4]
%                     * *                  *|*         ef2:[V4,V3]
%                    *   *                * | *        ef3:[V3,V5]
%                   *     *              *  |  *V5     ef4:[V5,V1]
%                 ec3      ec2          *   |RL/ *     ef5:[V1,V2]
%                 *   T=L   *          * L  | /  ef3   ef6:[V4,V5]
%                *           *        *     |/RR   *   ef7:[V5,V4]
%               V2*****ec1*****V3    V2*ef1*V4*ef2*V3 
        
NTc2 = length(splite2);  
if NTc2 > 0                                                     
    elem2edgec2f(splite2,1) = elem2edgef(splite2,2);
    elem2edgec2f(splite2,2) = elem2edgef(elemcl2r(splite2),1);
    elem2edgec2f(splite2,3) = elem2edgef(elemcl2r(splite2),2);
    elem2edgec2f(splite2,4) = elem2edgef(elemcl2r(elemcl2r(splite2)),3);
    elem2edgec2f(splite2,5) = elem2edgef(splite2,1);
    elem2edgec2f(splite2,6) = elem2edgef(splite2,3);
    elem2edgec2f(splite2,7) = elem2edgef(elemcl2r(splite2),3);

    elem2edgesignc2f(splite2,1) = elem2edgeSignf(splite2,2);
    elem2edgesignc2f(splite2,2) = elem2edgeSignf(elemcl2r(splite2),1);
    elem2edgesignc2f(splite2,3) = elem2edgeSignf(elemcl2r(splite2),2);
    elem2edgesignc2f(splite2,4) = elem2edgeSignf(elemcl2r(elemcl2r(splite2)),3);
    elem2edgesignc2f(splite2,5) = elem2edgeSignf(splite2,1);
    elem2edgesignc2f(splite2,6) = elem2edgeSignf(splite2,3);
    elem2edgesignc2f(splite2,7) = elem2edgeSignf(elemcl2r(splite2),3);

%              1     2    3     4     5      6      7
    locRij = [0.5  0.5   0      0     0      0    -0.5;...
              0     0    0.5   0.5    0      1    -0.5;...
              0     0    0      0     1     -1     0.5];
    for i = 1:3  
        for j = 1:7
            if locRij(i,j)~= 0
                Rij  = locRij(i,j)*elem2edgeSignc(splite2,i).*elem2edgesignc2f(splite2,j);
                ii(index+1:index+NTc2) = elem2edgec(splite2,i);
                jj(index+1:index+NTc2) = elem2edgec2f(splite2,j);
                ss(index+1:index+NTc2) = Rij/2;
                index = index + NTc2;
            end
        end
    end   
end

%**************************************************************************
%    CASE THREE: two bisections. left child is bisect.
%**************************************************************************
%
%                      V1                  V1          ef1:[V2,V4]
%                     * *                  *|*         ef2:[V4,V3]
%                    *   *              ef4 | *        ef3:[V3,V1]
%                   *     *              *LR|  *       ef4:[V1,V5]
%                 ec3      ec2         V5\  |  ef3     ef5:[V5,V2]
%                 *   T=L   *          *  \ |    *     ef6:[V1,V4]
%                *           *        * L  \| LR1 *    ef7:[V4,V5]
%               V2*****ec1*****V3    V2*ef1*V4*ef2*V3  

NTc3 = length(splite3);  
if NTc3 > 0    
    elem2edgec2f(splite3,1) = elem2edgef(elemcl2r(splite3),1);
    elem2edgec2f(splite3,2) = elem2edgef(elemcl2r1(splite3),3);
    elem2edgec2f(splite3,3) = elem2edgef(elemcl2r1(splite3),1);
    elem2edgec2f(splite3,4) = elem2edgef(splite3,2);
    elem2edgec2f(splite3,5) = elem2edgef(elemcl2r(splite3),3);
    elem2edgec2f(splite3,6) = elem2edgef(elemcl2r1(splite3),2);
    elem2edgec2f(splite3,7) = elem2edgef(elemcl2r(splite3),2);

    elem2edgesignc2f(splite3,1) = elem2edgeSignf(elemcl2r(splite3),1);
    elem2edgesignc2f(splite3,2) = elem2edgeSignf(elemcl2r1(splite3),3);
    elem2edgesignc2f(splite3,3) = elem2edgeSignf(elemcl2r1(splite3),1);
    elem2edgesignc2f(splite3,4) = elem2edgeSignf(splite3,2);
    elem2edgesignc2f(splite3,5) = elem2edgeSignf(elemcl2r(splite3),3);
    elem2edgesignc2f(splite3,6) = elem2edgeSignf(elemcl2r1(splite3),2);
    elem2edgesignc2f(splite3,7) = elem2edgeSignf(elemcl2r(splite3),2);

%              1     2    3    4      5      6      7
    locRij = [0.5   0.5   0     0     0     0     -0.5;...
              0     0     1     0     0     -1    0.5;...
              0     0     0    0.5    0.5    1   -0.5];

    for i = 1:3  
        for j = 1:7
            if locRij(i,j)~= 0
                Rij  = locRij(i,j)*elem2edgeSignc(splite3,i).*elem2edgesignc2f(splite3,j);
                ii(index+1:index+NTc3) = elem2edgec(splite3,i);
                jj(index+1:index+NTc3) = elem2edgec2f(splite3,j);
                ss(index+1:index+NTc3) = Rij/2;
                index = index + NTc3;
            end
        end
    end   
end

%**************************************************************************
%    CASE FOUR: three bisections. both children are bisected
%**************************************************************************
%                      V1                  V1          ef1:[V2,V4]
%                     * *                  *|*         ef2:[V4,V3]
%                    *   *                * | *        ef3:[V3,V6]
%                   *     *              *  |  *V6     ef4:[V6,V1]
%                 ec3      ec2         V5\ ef5 / *     ef5:[V1,V5]
%                 *   T=L   *          *  \ | /   *    ef6:[V5,V2]
%                *           *        *    \|/     *   ef7:[V4,V5]
%               V2*****ec1*****V3    V2*ef2*V4*ef3*V3  ef8:[V6,V4]  
%                                                      ef9:[V4,V1] 
NTc4 = length(splite4);  
if NTc4 >0
    elem2edgec2f(splite4,1) = elem2edgef(elemcl2r(splite4),1);
    elem2edgec2f(splite4,2) = elem2edgef(elemcl2r1(splite4),1);
    elem2edgec2f(splite4,3) = elem2edgef(elemcl2r1(splite4),2);
    elem2edgec2f(splite4,4) = elem2edgef(elemcl2r(elemcl2r1(splite4)),3);
    elem2edgec2f(splite4,5) = elem2edgef(splite4,2);
    elem2edgec2f(splite4,6) = elem2edgef(elemcl2r(splite4),3);
    elem2edgec2f(splite4,7) = elem2edgef(elemcl2r(splite4),2);
    elem2edgec2f(splite4,8) = elem2edgef(elemcl2r1(splite4),3);
    elem2edgec2f(splite4,9) = elem2edgef(splite4,1);

    elem2edgesignc2f(splite4,1) = elem2edgeSignf(elemcl2r(splite4),1);
    elem2edgesignc2f(splite4,2) = elem2edgeSignf(elemcl2r1(splite4),1);
    elem2edgesignc2f(splite4,3) = elem2edgeSignf(elemcl2r1(splite4),2);
    elem2edgesignc2f(splite4,4) = elem2edgeSignf(elemcl2r(elemcl2r1(splite4)),3);
    elem2edgesignc2f(splite4,5) = elem2edgeSignf(splite4,2);
    elem2edgesignc2f(splite4,6) = elem2edgeSignf(elemcl2r(splite4),3);
    elem2edgesignc2f(splite4,7) = elem2edgeSignf(elemcl2r(splite4),2);
    elem2edgesignc2f(splite4,8) = elem2edgeSignf(elemcl2r1(splite4),3);
    elem2edgesignc2f(splite4,9) = elem2edgeSignf(splite4,1);

%              1     2     3     4     5      6      7      8     9           
    locRij = [0.5   0.5    0     0     0      0    -0.5  -0.5   0;...
              0      0    0.5   0.5    0      0     0.5  -0.5   1;...
              0      0     0     0     0.5   0.5   -0.5   0.5  -1];

    for i = 1:3  
        for j = 1:9
            if locRij(i,j)~= 0
                Rij  = locRij(i,j)*elem2edgeSignc(splite4,i).*elem2edgesignc2f(splite4,j);
                ii(index+1:index+NTc4) = elem2edgec(splite4,i);
                jj(index+1:index+NTc4) = elem2edgec2f(splite4,j);
                ss(index+1:index+NTc4) = Rij/2;
                index = index + NTc4;
            end
        end
    end   
end

%% Modification for boundary edges
% remove zeros
idx = (ss == 0);
ss(idx) = [];
ii(idx) = [];
jj(idx) = [];
% Double the value for boundary edges
s = accumarray(elem2edgef(:), 1, [NEf 1]);
bdEdgeidxf = (s==1);
bdEdgeidxfinjj = bdEdgeidxf(jj);
ss(bdEdgeidxfinjj) = 2*ss(bdEdgeidxfinjj);
% Remove fixed boundary edges
if nargin>=5
    isBdEdgec = elem2edgec(bdFlagc(:) == 1);
    isBdEdgef = elem2edgef(bdFlagf(:) == 1);
    idx = ~(isBdEdgec(ii) | isBdEdgef(jj)); 
    Pro = sparse(jj(idx),ii(idx),ss(idx),NEf,NEc);
else
    Pro = sparse(jj,ii,ss,NEf,NEc);
end