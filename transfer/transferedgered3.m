function [Pro,freeEdgec] = transferedgered3(elemc,elemf,freeEdge)
%% TRANSFEREDGERED3 transfer operator of edge elements
%
%  Pro = transferedgered3(elemc,elemf) generates a prolongation matrix
%  between a tetrahedron mesh and its refinement by uniformrefine3.
%
% Written by Jie Zhou based on discussion with Long Chen on Aug 29 2013.
%
% See also: transferedgered
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if ~exist('freeEdge','var'), freeEdge = []; end 

%% Data structure on coarse and fine grids
[elem2edgec,edgec,elem2edgeSignc] = dof3edge(elemc);
[elem2edgef,edgef,elem2edgeSignf] = dof3edge(elemf);
NEc = size(edgec,1);  NEf = size(edgef,1); NTc = size(elemc,1);

%% Data structure from coarse to fine grids
% index map inside one element
elem2edgec2f     = zeros(NTc,25);
elem2edgec2f(:,1) = elem2edgef(1:NTc,1);   
elem2edgec2f(:,2) = elem2edgef(1:NTc,2); 
elem2edgec2f(:,3) = elem2edgef(1:NTc,3); 
elem2edgec2f(:,4) = elem2edgef(NTc+1:2*NTc,1); 
elem2edgec2f(:,5) = elem2edgef(NTc+1:2*NTc,4); 
elem2edgec2f(:,6) = elem2edgef(NTc+1:2*NTc,5); 
elem2edgec2f(:,7) = elem2edgef(2*NTc+1:3*NTc,2); 
elem2edgec2f(:,8) = elem2edgef(2*NTc+1:3*NTc,4); 
elem2edgec2f(:,9) = elem2edgef(2*NTc+1:3*NTc,6); 
elem2edgec2f(:,10) = elem2edgef(3*NTc+1:4*NTc,3); 
elem2edgec2f(:,11) = elem2edgef(3*NTc+1:4*NTc,5); 
elem2edgec2f(:,12) = elem2edgef(3*NTc+1:4*NTc,6); 
elem2edgec2f(:,13) = elem2edgef(1:NTc,4); 
elem2edgec2f(:,14) = elem2edgef(1:NTc,5); 
elem2edgec2f(:,15) = elem2edgef(NTc+1:2*NTc,2); 
elem2edgec2f(:,16) = elem2edgef(NTc+1:2*NTc,3);
elem2edgec2f(:,17) = elem2edgef(4*NTc+1:5*NTc,4);
elem2edgec2f(:,18) = elem2edgef(7*NTc+1:8*NTc,1);
elem2edgec2f(:,19) = elem2edgef(7*NTc+1:8*NTc,2);
elem2edgec2f(:,20) = elem2edgef(7*NTc+1:8*NTc,3);
elem2edgec2f(:,21) = elem2edgef(6*NTc+1:7*NTc,4); 
elem2edgec2f(:,22) = elem2edgef(6*NTc+1:7*NTc,5);
elem2edgec2f(:,23) = elem2edgef(7*NTc+1:8*NTc,4);
elem2edgec2f(:,24) = elem2edgef(7*NTc+1:8*NTc,5);
elem2edgec2f(:,25) = elem2edgef(7*NTc+1:8*NTc,6); 

% sign map inside one element
elem2edgeSignc2f = sparse(NTc,25);
elem2edgeSignc2f(:,1) = elem2edgeSignf(1:NTc,1);   
elem2edgeSignc2f(:,2) = elem2edgeSignf(1:NTc,2); 
elem2edgeSignc2f(:,3) = elem2edgeSignf(1:NTc,3); 
elem2edgeSignc2f(:,4) = elem2edgeSignf(NTc+1:2*NTc,1); 
elem2edgeSignc2f(:,5) = elem2edgeSignf(NTc+1:2*NTc,4); 
elem2edgeSignc2f(:,6) = elem2edgeSignf(NTc+1:2*NTc,5); 
elem2edgeSignc2f(:,7) = elem2edgeSignf(2*NTc+1:3*NTc,2); 
elem2edgeSignc2f(:,8) = elem2edgeSignf(2*NTc+1:3*NTc,4); 
elem2edgeSignc2f(:,9) = elem2edgeSignf(2*NTc+1:3*NTc,6); 
elem2edgeSignc2f(:,10) = elem2edgeSignf(3*NTc+1:4*NTc,3); 
elem2edgeSignc2f(:,11) = elem2edgeSignf(3*NTc+1:4*NTc,5); 
elem2edgeSignc2f(:,12) = elem2edgeSignf(3*NTc+1:4*NTc,6); 
elem2edgeSignc2f(:,13) = elem2edgeSignf(1:NTc,4); 
elem2edgeSignc2f(:,14) = elem2edgeSignf(1:NTc,5); 
elem2edgeSignc2f(:,15) = elem2edgeSignf(NTc+1:2*NTc,2); 
elem2edgeSignc2f(:,16) = elem2edgeSignf(NTc+1:2*NTc,3);
elem2edgeSignc2f(:,17) = elem2edgeSignf(4*NTc+1:5*NTc,4);
elem2edgeSignc2f(:,18) = elem2edgeSignf(7*NTc+1:8*NTc,1);
elem2edgeSignc2f(:,19) = elem2edgeSignf(7*NTc+1:8*NTc,2);
elem2edgeSignc2f(:,20) = elem2edgeSignf(7*NTc+1:8*NTc,3);
elem2edgeSignc2f(:,21) = elem2edgeSignf(6*NTc+1:7*NTc,4); 
elem2edgeSignc2f(:,22) = elem2edgeSignf(6*NTc+1:7*NTc,5);
elem2edgeSignc2f(:,23) = elem2edgeSignf(7*NTc+1:8*NTc,4);
elem2edgeSignc2f(:,24) = elem2edgeSignf(7*NTc+1:8*NTc,5);
elem2edgeSignc2f(:,25) = elem2edgeSignf(7*NTc+1:8*NTc,6); 

%% Local restriction operator
c1 = 0.5;
c2 = 0.25;
locRij = zeros(6,25);
locRij(6,25) = -c2; locRij(5,25) =  c2; locRij(4,25) = c2;
locRij(6,24) =  c2; locRij(5,24) =  c2; locRij(4,24) = c2;
locRij(6,23) =  c2; locRij(5,23) =  c2; locRij(4,23) = -c2;
locRij(6,22) = -c2; locRij(3,22) =  c2; locRij(2,22) = c2;
locRij(5,21) = -c2; locRij(3,21) =  c2; locRij(1,21) = c2;
locRij(6,20) =  c2; locRij(3,20) =  c2; locRij(2,20) = c2;
locRij(6,19) =  c2; locRij(4,19) = -c2; locRij(3,19) = c2;  locRij(1,19)= c2;
locRij(4,18) = -c2;locRij(2,18) =  c2;locRij(1,18) = c2;
locRij(6,17) =  c2;locRij(3,17) =  c2;locRij(2,17) = -c2;
locRij(5,16) =  c2;locRij(3,16) =  c2;locRij(1,16) = c2;
locRij(4,15) =  c2;locRij(2,15) =  c2;locRij(1,15) = c2;
locRij(5,14) =  c2;locRij(3,14) =  c2;locRij(1,14) = -c2;
locRij(4,13) =  c2;locRij(2,13) =  c2;locRij(1,13) = -c2;
locRij(6,12) =  c1;
locRij(5,11) =  c1;
locRij(3,10) =  c1;
locRij(6,9)  =  c1;
locRij(4,8)  =  c1;
locRij(2,7)  =  c1;
locRij(5,6)  =  c1;
locRij(4,5)  =  c1;
locRij(1,4)  =  c1;
locRij(3,3)  =  c1;
locRij(2,2)  =  c1;
locRij(1,1)  =  c1;

%% Assemble the prolongation matrix
ii = zeros(25*NTc,1); jj = zeros(25*NTc,1); ss = zeros(25*NTc,1);multiplity = zeros(25*NTc,1);
index = 0;
for i = 1:6  
    for j = 1:25
        if locRij(i,j)~=0  
            Rij = locRij(i,j)*elem2edgeSignc(:,i).*elem2edgeSignc2f(:,j);
            ii(index+1:index+NTc) =  elem2edgec(:,i);
            jj(index+1:index+NTc) =  elem2edgec2f(:,j);
            ss(index+1:index+NTc) =  Rij;
            multiplity(index+1:index+NTc) = ones(NTc,1);
            index = index + NTc;
         end
    end
end

%% Modification for duplicated edges
% We work elementwise and thus an edge is repeated several times.
Pro     =  sparse(jj,ii,ss,NEf,NEc);
Multi   =  sparse(jj,ii,multiplity,NEf,NEc);
[i1,j1,s1] = find(Pro);
[tempi1,tempi2,s2]   = find(Multi); %#ok<*ASGLU>
Pro = sparse(i1,j1,s1./s2,NEf,NEc);

%% Truced to free edges
if ~isempty(freeEdge)
    isFixedEdge  = true(NEf,1);
    isFixedEdge(freeEdge) = false;
    isFixedEdgec = false(NEc,1);
    fixedEdgeNode  = edgef(isFixedEdge,:);
    index = fixedEdgeNode(:)-Nc;
    index(index<=0) = [];
    isFixedEdgec(index) = true;
    freeEdgec = find(~isFixedEdgec);
    Pro = Pro(freeEdge,freeEdgec);
end