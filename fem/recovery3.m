function RDu = recovery3(node,elem,Du,volume)
%% RECOVERY3 recovery a piecewise constant function to a piecewise linear one in 3-D.
%
% RDu = recovery3(node,elem,Du,volume) compute a P1 approximation u using area
% weighted average of the piecewise constant gradient Du usually given by
% [Du,area] = gradu3(node,elem,u).
%
% See also gradu3, estimaterecovery, recovery3
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

N = size(node,1);
dudxVolume = volume.*Du(:,1);
dudyVolume = volume.*Du(:,2);
dudzVolume = volume.*Du(:,3);
patchVolume = accumarray(elem(:),repmat(volume,4,1), [N 1]); 
dudxVolume = accumarray(elem(:),repmat(dudxVolume,4,1),[N 1]);
dudyVolume = accumarray(elem(:),repmat(dudyVolume,4,1),[N 1]);
dudzVolume = accumarray(elem(:),repmat(dudzVolume,4,1),[N 1]);
dudx = dudxVolume./patchVolume;
dudy = dudyVolume./patchVolume;
dudz = dudzVolume./patchVolume;
RDu = [dudx, dudy, dudz];