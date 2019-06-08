function RDu = recovery(node,elem,Du,area)
%% RECOVERY recovery a piecewise constant function to a piecewise linear one.
%
% RDu = recovery(node,elem,Du,area) compute a P1 approximation u using area
% weighted average of the piecewise constant gradient Du usually given by
% [Du,area] = gradu(node,elem,u).
%
% See also gradu, estimaterecovery, recovery3
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

N = size(node,1);
dudxArea = area.*Du(:,1);
dudyArea = area.*Du(:,2);
patchArea = accumarray(elem(:),[area;area;area], [N 1]); 
dudxArea = accumarray(elem(:),[dudxArea;dudxArea;dudxArea],[N 1]);
dudyArea = accumarray(elem(:),[dudyArea;dudyArea;dudyArea],[N 1]);
dudx = dudxArea./patchArea;
dudy = dudyArea./patchArea;
RDu = [dudx, dudy];
