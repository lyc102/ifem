%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function d = dCircle(P,xc,yc,r)
d=sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
d=[d,d];
%-------------------------------------------------------------------------%