function markedElem = mark(elem,eta,theta,method)
% MARK mark element.
%
% markedElem = mark(elem,eta,theta) mark a subset of elements by Dorfler
% marking strategy. It returns an array of indices of marked elements
% markedElem such that sum(eta(markedElem)^2) > theta*sum(eta^2).
%
% markedElem = mark(elem,eta,theta,'max') choose markedElem such that
% eta(markedElem) > theta*max(eta).
%
% markedElem = mark(elem,eta,theta,'COARSEN') choose markedElem such that
% eta(markedElem) < theta*max(eta).
%
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.

NT = size(elem,1); isMark = false(NT,1);
if ~exist('method','var'), method = 'L2'; end  % default marking is L2 based
switch upper(method)
    case 'MAX'
        isMark(eta>theta*max(eta))=1;
    case 'COARSEN'
        isMark(eta<theta*max(eta))=1;
    case 'L2'
        [sortedEta,idx] = sort(eta.^2,'descend'); 
        x = cumsum(sortedEta);
        isMark(idx(x < theta* x(NT))) = 1;
        isMark(idx(1)) = 1;
end
markedElem = uint32(find(isMark==true));