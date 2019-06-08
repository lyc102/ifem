function [Pro,freeEdgec] = transferedgered(elemc,elemf,freeEdge)
%% TRANSFEREDGERED transfer operator of edge elements
%
% Pro = transferedgered(elemc,elemf) generates a prolongation matrix
%  between a triangulationand its refinement by uniformrefine.
%
% Authors: Long Chen, Ming Wang, Jie Zhou. 
%
% See also: transferedgered3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if ~exist('freeEdge','var'), freeEdge = []; end 

%% Data structure on coarse and fine grids
[elem2edgec,edgec,elem2edgeSignc] = dofedge(elemc);
[elem2edgef,edgef,elem2edgeSignf] = dofedge(elemf);
elem2edgeSignf = double(elem2edgeSignf);
elem2edgeSignc = double(elem2edgeSignc);
NEc = size(edgec,1);  NEf = size(edgef,1); NTc = size(elemc,1);  

%% Data structure from coarse to fine grids
% index map inside one element
elem2edgec2f = zeros(NTc,9);
elem2edgec2f(:,1) = elem2edgef(NTc+1:2*NTc,1);
elem2edgec2f(:,2) = elem2edgef(2*NTc+1:3*NTc,1);
elem2edgec2f(:,3) = elem2edgef(2*NTc+1:3*NTc,2);
elem2edgec2f(:,4) = elem2edgef(1:NTc,2);
elem2edgec2f(:,5) = elem2edgef(1:NTc,3);
elem2edgec2f(:,6) = elem2edgef(NTc+1:2*NTc,3);
elem2edgec2f(:,7) = elem2edgef(1:NTc,1);
elem2edgec2f(:,8) = elem2edgef(NTc+1:2*NTc,2);
elem2edgec2f(:,9) = elem2edgef(2*NTc+1:3*NTc,3);
% sign map inside one element
elem2edgeSignc2f = zeros(NTc,9);
elem2edgeSignc2f(:,1) = elem2edgeSignf(NTc+1:2*NTc,1);
elem2edgeSignc2f(:,2) = elem2edgeSignf(2*NTc+1:3*NTc,1);
elem2edgeSignc2f(:,3) = elem2edgeSignf(2*NTc+1:3*NTc,2);
elem2edgeSignc2f(:,4) = elem2edgeSignf(1:NTc,2);
elem2edgeSignc2f(:,5) = elem2edgeSignf(1:NTc,3);
elem2edgeSignc2f(:,6) = elem2edgeSignf(NTc+1:2*NTc,3);
elem2edgeSignc2f(:,7) = elem2edgeSignf(1:NTc,1);
elem2edgeSignc2f(:,8) = elem2edgeSignf(NTc+1:2*NTc,2);
elem2edgeSignc2f(:,9) = elem2edgeSignf(2*NTc+1:3*NTc,3);

% Local edge labling
%
%                V1
%                *
%               * *
%              *   *
%             5     4
%            *       *
%           *         *
%          *           *
%         E3------7------E2
%        * \             /*
%       *   \           /  *
%      *     \         /    *
%     6       8       9      3
%    *         \     /        *
%   *           \   /          *
%  *             \ /            *
% V2 **** 1 **** E1 ***** 2 **** V3

%% Local restriction operator
%          1     2     3     4      5     6      7      8     9           
locRij = [0.5   0.5    0     0      0     0     0.25   -0.25   -0.25;...
          0      0   0.5     0.5    0     0    -0.25    0.25   -0.25;...
          0      0     0     0     0.5   0.5   -0.25   -0.25    0.25];

%% Assemble the prolongation matrix
ii = zeros(15*NTc,1); jj = zeros(15*NTc,1); ss = zeros(15*NTc,1);
index = 0;
for i = 1:3  
    for j = 1:9
        if locRij(i,j)~=0
            Rij = locRij(i,j)*elem2edgeSignc(:,i).*elem2edgeSignc2f(:,j);
            ii(index+1:index+NTc) = elem2edgec(:,i);
            jj(index+1:index+NTc) = elem2edgec2f(:,j);
            if j<7
                ss(index+1:index+NTc) = 0.5*Rij;
            % Every interior edges of the coarse grid will be computed
            % twice. When assembling, the value is halfed. Values for
            % boundary edges will be doubled later.
            else
                ss(index+1:index+NTc) = Rij;                
            end
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
    isFixedEdgec = false(NEc,1);
    fixedEdgeNode  = edgef(isFixedEdge,:);
    Nc = max(elemc(:));
    index = fixedEdgeNode(:)- Nc; % nodes in the fine grid only is on edges
    index(index<=0) = [];
    isFixedEdgec(index) = true;
    freeEdgec = find(~isFixedEdgec);
    Pro = Pro(freeEdge,freeEdgec);
end