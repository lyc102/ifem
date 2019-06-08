function elemType = updatetype(elemType,brother)
%% UPDATETYPE update the type of triangles.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if (size(brother,1)==0), return; end
oldNT = length(elemType); 
newNT = max(brother(:,2)); % brother(t,1) < brother(t,2)
if oldNT < newNT % coarse grid to fine grid
    elemType(newNT) = 0;    
    firstBrother = zeros(oldNT,1,'uint32');
    idx = (brother(:,1) <= oldNT);
    t = brother(idx,1);
    tt = brother(idx,2);
    firstBrother(t(end:-1:1)) = tt(end:-1:1);
    % update first bisected elements from 1:oldNT first
    elemType(firstBrother(t)) = ~elemType(t);
    elemType(t) = ~elemType(t);
    % updated second bisected elements
    idx = (brother(:,2) > (oldNT + length(find(firstBrother))));
    elemType(brother(idx,2)) = ~elemType(brother(idx,1));
    elemType(brother(idx,1)) = elemType(brother(idx,2));
else % fine grid to coarse grid
    elemType(brother(:,1)) = ~elemType(brother(:,1));
    elemType(brother(:,2)) = [];
end