function [p,t] = cleanup(p,t,fd,pfix,varargin)
%% CLEANUP delete isolated vertices.
%
% [p,t] = cleanup(p,t) rmove isolate points in the triangulation
% (p,t). Isolate points are interiori vertices with valence 3 and
% boundary vertices with valence 2. Here the valence of a vertex is defined
% as the triangles containing that vertex.
%
% Example
%   load airfoilperturb
%   [node,elem] = cleanup(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

N = size(p,1); NT = size(t,1); Nfix = size(pfix,1);

bdNode = findboundary(t);
isBdNode = false(N,1);
isBdNode(bdNode) = true;
intNode = find(~isBdNode);

valence = accumarray(t(:),ones(3*NT,1),[N 1]);

isIsoNode = false(N,1);
isIsoNode(bdNode(valence(isBdNode)==2)) = true;
isIsoNode(intNode(valence(~isBdNode)<=4)) = true;
isIsoNode(1:Nfix) = false;
isoNode = find(isIsoNode);

isIsoElem = isIsoNode(t(:,1)) | isIsoNode(t(:,2)) | isIsoNode(t(:,3));
tt = t(isIsoElem,:);
q = simpqual(p, tt);
Ntt = size(tt,1);
t2v = sparse([1:Ntt,1:Ntt,1:Ntt], tt(1:Ntt,:), 1, Ntt, N);
for i = 1:length(isoNode)
    ring = find(t2v(:,isoNode(i)));
    qi = min(q(ring)); %#ok<FNDSB>
    if qi>0.65
        isIsoNode(isoNode(i)) = false;
    end
end

isNotUsed = true(N,1);
isNotUsed(t(:)) = false;
rmpt = isIsoNode | isNotUsed;

if find(rmpt)
%     findnode(p,isIsoNode,'noindex','color','r')
    p(rmpt,:) = [];
%     fprintf('%d points are deleted \n', sum(rmpt))
    geps = 0.001/sqrt(NT); 
    t = delaunayn(p); 
    center = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3; 
    t = t(feval(fd,center,varargin{:})<-geps,:);
    t = fixorder(p,t);
end