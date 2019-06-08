function [interfaceElem,coarsenElem,isInterfaceElem] = markinterface(elem,phiValue)
%% MARKINTERFACE mark interface elements
%
% Interface element is defined as the element whose vertices are not in the
% same side of the interface. It is found by checking the sign of the
% values of the level set function phi.
%
% See also: interfacemesh
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

s1 = msign(phiValue(elem(:,1)));
s2 = msign(phiValue(elem(:,2)));
s3 = msign(phiValue(elem(:,3)));
eta1 = abs(s1+s2+s3);
eta2 = abs(s1)+abs(s2)+abs(s3);
isInterfaceElem = (eta1 == 1 & eta2 == 3)|(eta1 == 0 & eta2 == 2) ;
% case: eta1 == 1 & eta2 == 3
%  eta2 = 3: no interface nodes
%  eta1 = 1: at least one edge cross the interface
% case: eta1 == 0 & eta2 == 2
%  eta2 = 2: one vertex is on the interface
%  eta1 = 0: at least one edge cross the interface
% Note that if eta2 <2, then there are at least two interface nodes and the
% element is not an interface element
interfaceElem  = find(isInterfaceElem); % elements cross interface 
coarsenElem = find(~isInterfaceElem);   % elements away from interface