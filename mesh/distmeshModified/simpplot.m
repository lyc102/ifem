function simpplot(p,t,expr,bcol,icol)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.

if nargin<4, bcol=.9*ones(1,3); end
if nargin<5, icol=.6*ones(1,3); end

dim=size(p,2);
switch dim
 case 2
  trimesh(t,p(:,1),p(:,2),0*p(:,1),'facecolor','none','edgecolor','k');
  view(2)
  axis equal
  axis off
 case 3
  tri1=surftri(p,t);
  if nargin>2 & ~isempty(expr)
    incl=find(eval(expr));
    t=t(any(ismember(t,incl),2),:);
    tri1=tri1(any(ismember(tri1,incl),2),:);
    tri2=surftri(p,t);
    tri2=setdiff(tri2,tri1,'rows');
    h=trimesh(tri2,p(:,1),p(:,2),p(:,3));
    set(h,'facecolor',icol,'edgecolor','k');
    hold on
  end
  h=trimesh(tri1,p(:,1),p(:,2),p(:,3));
  hold off
  set(h,'facecolor',bcol,'edgecolor','k');
  axis equal
  cameramenu
 otherwise
  error('Unimplemented dimension.');
end
