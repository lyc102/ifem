%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function d = dIntersect(d1,d2) % max(d1,d2)
d=[d1(:,1:(end-1)),d2(:,1:(end-1))];
d=[d,max(d1(:,end),d2(:,end))];
%-------------------------------------------------------------------------%