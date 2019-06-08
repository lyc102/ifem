function [p,t]=fixmesh(p,t)
%FIXMESH  Remove duplicated/unused nodes and fix element orientation.
%   [P,T]=FIXMESH(P,T)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.

snap=max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[foo,ix,jx]=unique(round(p/snap)*snap,'rows');
p=p(ix,:);
t=jx(t);

[pix,ix,jx]=unique(t);
t=reshape(jx,size(t));
p=p(pix,:);

flip=simpvol(p,t)<0;
t(flip,[1,2])=t(flip,[2,1]);
