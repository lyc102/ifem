function [elem,HB,bdFlag] = uniformcoarsenred(elem,bdFlag)
%% UNIFORMCOARSENRED uniform coarsening of red refinement
%
% [elem,HB] = uniformcoarsenred(elem) remove grid points added by uniform
% refinement. See the illustration below:
%
% 3
% | \
% 5- 4    t3
% |\ |\    t4
% 1-6 -2  t1 t2
%
% It is mainly used to get multilevel decomposition in multigrid methods.
% The input matrix elem stands for the fine mesh and the output one for the
% coarse mesh. The HB records the hierarchical structure of added points
% going from the coarse to the fine mesh such that HB(:,2:3) are two parent
% nodes of HB(:,1).
%
%   See also: uniformcoarsen, coarsen, bisect, uniformcoarsen3, mg
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

HB = [];
NT = size(elem,1);
if mod(NT,4)==0
    NTc = NT/4; % number of triangles in the coarse grid
else
%     display('Not from red refinement');
    return
end

%% Find points
t1 = 1:NTc; t2 = t1+NTc; t3 = t2+NTc; t4 = t3+NTc;
if any(elem(t1,2)~=elem(t2,1)) || any(elem(t1,3)~=elem(t3,1)) || ...
   any(elem(t4,1)~=elem(t2,3)) || any(elem(t4,2)~=elem(t3,1))
%     display('Not from red refinement');
    return
end
p1 = elem(t1,1);
p2 = elem(t2,2);
p3 = elem(t3,3);
p4 = elem(t4,1);
p5 = elem(t1,3);
p6 = elem(t1,2);

%% Update and remove triangles
elem(t1,:) = [p1 p2 p3];
elem = elem(t1,:);

%% Record HB
HB(p6,:) = [p6 p1 p2];
HB(p4,:) = [p4 p2 p3];
HB(p5,:) = [p5 p1 p3];
Nc = max(elem(:));
HB = HB(Nc+1:end,:);

%% Update boundary edges
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlag(t1,1) = bdFlag(t2,1);
    bdFlag = bdFlag(t1,:);
else
    bdFlag = [];
end