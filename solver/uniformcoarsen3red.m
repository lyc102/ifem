function [elem,HBmesh,HB,bdFlag] = uniformcoarsen3red(elem,HBmesh,bdFlag)
%% UNIFORMCOARSEN3RED uniform coarsening of red refinement
%
% [elem,HB] = uniformcoarsen3red(elem) remove grid points added by uniform
% refinement in 3-D. 
%
% It is mainly used to get multilevel decomposition in multigrid methods.
% The input matrix elem stands for the fine mesh and the output one for the
% coarse mesh. The HB records the hierarchical structure of added points
% going from the coarse to the fine mesh such that HB(:,2:3) are two parent
% nodes of HB(:,1).
%
%   See also: uniformrefine3, mg, uniformcoarsen, uniformcoarsen3
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

HB = [];
NT = size(elem,1);
if mod(NT,8)==0
    NTc = NT/8; % number of triangles in the coarse grid
else
%     display('Not from red refinement');
    return
end

%% Find points
t1 = 1:NTc; t2 = t1+NTc; t3 = t2+NTc; t4 = t3+NTc;
if any(elem(t1,2)~=elem(t2,1)) || any(elem(t1,3)~=elem(t3,1)) || ...
   any(elem(t4,1)~=elem(t1,4)) || any(elem(t4,2)~=elem(t2,4))
%     display('Not from red refinement');
    return
end
p1 = elem(t1,1);
p2 = elem(t2,2);
p3 = elem(t3,3);
p4 = elem(t4,4);
p5 = elem(t1,2);
p6 = elem(t1,3);
p7 = elem(t1,4);
p8 = elem(t2,3);
p9 = elem(t2,4);
p10 = elem(t3,4);

%% Remove tetrahedron
elem(t1,:) = [p1 p2 p3 p4];
elem = elem(t1,:);

%% Record HB
HB(p5,:) = [p5 p1 p2];
HB(p6,:) = [p6 p1 p3];
HB(p7,:) = [p7 p1 p4];
HB(p8,:) = [p8 p2 p3];
HB(p9,:) = [p9 p2 p4];
HB(p10,:) = [p10 p3 p4];
Nc = max(elem(:));
HB = HB(Nc+1:end,:);

%% Update boundary edges
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlag(t1,1) = bdFlag(t2,1);
    bdFlag = bdFlag(t1,:);
else
    bdFlag = [];
end

%% Update HBmesh
if exist('HBmesh','var') && ~isempty(HBmesh)
    HBmesh = HBmesh(1:Nc,:);
else
    HBmesh = [];
end