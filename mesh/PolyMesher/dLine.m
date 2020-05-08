%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function d = dLine(P,x1,y1,x2,y2)
% By convention, a point located at the left hand side of the line
% is inside the region and it is assigned a negative distance value.
a = [x2-x1,y2-y1]; a = a/norm(a);
b = [P(:,1)-x1,P(:,2)-y1];
d = b(:,1)*a(2) - b(:,2)*a(1);
d = [d,d];
%-------------------------------------------------------------------------%