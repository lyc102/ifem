function [face, face2elem] = getpolymesh(elem, interfaceData)

elemIdx = interfaceData.elemIdx;
elem = elem(elemIdx ~= 0, :);
idx = find(elemIdx ~= 0);

face =  [elem(:,1) elem(:,4) elem(:,3) elem(:,2);...
         elem(:,6) elem(:,7) elem(:,8) elem(:,5); ...
         elem(:,2) elem(:,3) elem(:,7) elem(:,6);...
         elem(:,1) elem(:,5) elem(:,8) elem(:,4);...
         elem(:,1) elem(:,2) elem(:,6) elem(:,5);...
         elem(:,4) elem(:,8) elem(:,7) elem(:,3);...
         interfaceData.face];
face2elem = [idx;idx;idx;idx;idx;idx;interfaceData.face2elem];

NP = max(face2elem);
isPoly = false(NP, 1);
isPoly(face2elem) = true;
NNP = sum(isPoly);

idx2idx = zeros(NP, 1);
idx2idx(isPoly) = 1:NNP;
face2elem = idx2idx(face2elem);

end



