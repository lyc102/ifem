function [Pro,freeEdgec] = transferedgecoarsen(elemc,elemf,tree,freeEdge)
%% TRANSFEREDGECOARSEN transfer operator between bisected grids
%
%  Works for biscetion grids.
%
%  input elemc,elemf,brother
%  ouput Pro,Res
%   
% Added by Jie Zhou. Rewritten by Long Chen.

if ~exist('freeEdge','var'),freeEdge = []; end

%% Data structure
% elem2edge and elem2edgeSign
[elem2edgec,edgec,elem2edgeSignc] = dofedge(elemc);
[elem2edgef,edgef,elem2edgeSignf] = dofedge(elemf);
% number of edges and elements
NEc = size(edgec,1);  NEf = size(edgef,1); 
NTc = size(elemc,1);  NTf = size(elemf,1);
% tree 
%   t
%  / \
% tl tr
t = tree(:,1);  % parent element in the coarse grid 
tl = tree(:,2); % left child
tr = tree(:,3); % right child
% index map
inCoarse = true(NTf,1);
inCoarse(tr) = false;      % tr is removed from the coarsening
elemf2c = zeros(NTf,1);    
elemf2c(inCoarse) = 1:NTc; % fine to coarse index map of elements 
% elemc2f = find(elemf2c); % coarse to fine index map of elements

%% Assemble the prolongation matrix
% case 1: refined once
%                      V1                  V1       ef1:[V1,V2]
%                     * *                  *|*      ef2:[V2,V4]
%                    *   *                * | *     ef3:[V4,V3]
%                   *     *              *  |  *    ef4:[V3,V1]
%                 ec3      ec2         ef1 ef5 ef4  ef5:[V1,V4] 
%                 *    T    *          * L  | R  *  
%                *           *        *     |      *
%              V2 *** ec1 *** V3    V2*ef2* V4*ef3*V3  
% preallocation
NTc1 = length(t);  
ii = zeros(6*NTc1,1); jj = zeros(6*NTc1,1); ss = zeros(6*NTc1,1);
index = 0;
if NTc1 > 0
    elem2edgec2f(t,1) = elem2edgef(tl,1);
    elem2edgec2f(t,2) = elem2edgef(tl,2);
    elem2edgec2f(t,3) = elem2edgef(tr,3);
    elem2edgec2f(t,4) = elem2edgef(tr,1);
    elem2edgec2f(t,5) = elem2edgef(tr,2);
    elem2edgeSignc2f(t,1) = elem2edgeSignf(tl,1);
    elem2edgeSignc2f(t,2) = elem2edgeSignf(tl,2);
    elem2edgeSignc2f(t,3) = elem2edgeSignf(tr,3);
    elem2edgeSignc2f(t,4) = elem2edgeSignf(tr,1);
    elem2edgeSignc2f(t,5) = elem2edgeSignf(tr,2);
%             1     2     3    4     5
    locRij = [0   0.5   0.5    0      0  ;... % 1
              0     0     0    1   -0.5  ;... % 2
              1     0     0    0    0.5 ];    % 3
    [i,j,locRij] = find(locRij); % change to one dimensional array
    for k = 1:length(i);
        Rij  = locRij(k)*elem2edgeSignc(t,i(k)).*elem2edgeSignc2f(t,j(k));
        ii(index+1:index+NTc1) = elem2edgec(t,i(k));
        jj(index+1:index+NTc1) = elem2edgec2f(t,j(k));
        if j(k) == 5
            ss(index+1:index+NTc1) = Rij;            
        else
            ss(index+1:index+NTc1) = Rij/2;
        end
        index = index + NTc1;
    end
end
% case 0: non-refined
isNotRefined = true(NTf,1);
isNotRefined(tl) = false;  % the left one is the old index but refined
isNotRefined(tr) = false;  % new added elem is of course refined
notRefinedidx = find(isNotRefined); % not-refined element in the fine mesh
notRefinedCoarseIdx = elemf2c(notRefinedidx); % idx in the coarse mesh
NTc0 = length(notRefinedCoarseIdx);
if NTc0 > 0
    i0 = elem2edgec(notRefinedCoarseIdx,:);
    j0 = elem2edgef(notRefinedidx,:); 
    ii = [ii; i0(:)];
    jj = [jj; j0(:)];
    ss = [ss; ones(3*NTc0,1)/2];
end

%% Modification for boundary edges
% Double the value for boundary edges
s = accumarray(elem2edgef(:), 1, [NEf 1]);
bdEdgeidxf = (s==1);
bdEdgeidxfinjj = bdEdgeidxf(jj);
ss(bdEdgeidxfinjj) = 2*ss(bdEdgeidxfinjj);
Pro = sparse(double(jj),double(ii),ss,NEf,NEc);
% Remove fixed boundary edges
if ~isempty(freeEdge)
    isFixedEdgef = true(NEf,1);
    isFixedEdgef(freeEdge) = false;
    [i,j] = find(Pro(isFixedEdgef,:)); %#ok<ASGLU>
    isFreeEdgec = true(NEc,1);
    isFreeEdgec(j) = false;
    freeEdgec = find(isFreeEdgec);
    Pro = Pro(freeEdge,freeEdgec);
else
    freeEdgec = [];
end