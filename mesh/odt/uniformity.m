function u=uniformity(p,t,fh,varargin)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
if size(p,2) == 3
    trep = TriRep(t,p(:,1),p(:,2),p(:,3));
elseif size(p,2) == 2
    trep = TriRep(t,p(:,1),p(:,2));
end
[pc,r] = circumcenters(trep);
if isnumeric(fh)
    hc = fh;
else
    hc = feval(fh,pc,varargin{:});
end
sz=r./hc;
u=std(sz)/mean(sz);
