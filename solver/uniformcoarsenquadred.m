function [elem,HB] = uniformcoarsenquadred(elem)
%% UNIFORMCOARSEQUADNRED uniform coarsening of red refinement
%
% [elem,HB] = uniformcoarsenred(elem) remove grid points added by uniform
% refinement. See the illustration below:
%
% 4  - 7 -  3
% | t4 | t3 |
% 8 -  9 -  6
% | t1 | t2 |
% 1 -  5 -  2
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
   any(elem(t4,1)~=elem(t1,4)) || any(elem(t4,2)~=elem(t3,1))
%     display('Not from red refinement');
    return
end
p1 = elem(t1,1);
p2 = elem(t2,2);
p3 = elem(t3,3);
p4 = elem(t4,4);
p5 = elem(t1,2);
p6 = elem(t2,3);
p7 = elem(t3,4);
p8 = elem(t4,1);
p9 = elem(t1,3);

%% Remove quad
elem(t1,:) = [p1 p2 p3 p4];
elem = elem(t1,:);

%% Record HB
HB(p9,:) = [p9 p1 p3];
HB(p5,:) = [p5 p1 p2];
HB(p6,:) = [p6 p2 p3];
HB(p7,:) = [p7 p4 p3];
HB(p8,:) = [p8 p1 p4];
Nc = max(elem(:));
HB = HB(Nc+1:end,:);