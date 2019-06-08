function [node,elem,indexMap] = rmisopoint(node,elem,gflag)
%% RMISOPOINT remove isolate points 
%
% [node,elem] = rmisopoint(node,elem) removes isolate points in the
% triangulation (node,elem). Isolate points are interiori vertices with
% valence 3 and boundary vertices with valence 2. Here the valence of a
% vertex is defined as the triangles containing that vertex.
%
% [node,elem,indexMap] = rmisopoint(node,elem) also return the index map
% between the input and output vertices. Since some vertices could be
% deleted in the procedure, the indexMap is important for the interpolation
% of functions defined on these vertices. See nodeinterpolate.
%
% [node,elem] = rmisopoint(node,elem,1) will plot isolated points in red.
%
% Example
%   load airfoilperturb
%   [node,elem] = rmisopoint(node,elem);
%
% See also  bdsmoothing, edgeswap, nodeinterpolate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


N = size(node,1); NT = size(elem,1);
%% Find boundary nodes
bdNode = findboundary(elem);
isBdNode = false(N,1);
isBdNode(bdNode) = true;
isIntNode = ~isBdNode;
intNode = find(~isBdNode);

%% Compute valence and find isolate points
valence = accumarray(elem(:),ones(3*NT,1),[N 1]);
isIsoNode = false(N,1);
isIsoNode(bdNode(valence(isBdNode)==2)) = true;
isIsoNode(intNode(valence(isIntNode)==3)) = true;
isoNode = find(isIsoNode);

%% Construct node star of isolated points
isIsoElem = isIsoNode(elem(:,1))|isIsoNode(elem(:,2))|isIsoNode(elem(:,3));
isoElem = find(isIsoElem);
tt = elem(isIsoElem,:);
q = simpqual(node,tt);
Ntt = size(tt,1);
tt2v = sparse([1:Ntt,1:Ntt,1:Ntt], tt(1:Ntt,:), 1, Ntt, N);

%% Remove isolate points
for i = 1:length(isoNode)
    pi = isoNode(i);
    ring = find(tt2v(:,pi));
    qi = min(q(ring));
    if qi>0.65  % don't remove this point
        isIsoNode(isoNode(i)) = false;
    else % remove this isolate point
        allpt = tt(ring,:);
        linkpt = setdiff(unique(allpt(:)),pi);
        if isIntNode(pi) % interior isolate points
            newt = fixorder(node,linkpt');
            tt(ring(1),:) = newt;
            tt(ring(2:end),:) = 0;
            % isoElem is used as an index mapping from tt to t
            elem(isoElem(ring(1)),:) = newt;
            elem(isoElem(ring(2:end)),1) = 0;
        elseif (isBdNode(pi) && isempty(find(linkpt == 0,1)))
        % Corner/feature points could be isolate nodes. Remove these
        % points will change the shape (and thus the area) of the mesh. 
            newt = fixorder(node,linkpt');
            newarea = simplexvolume(node,newt);
            oldarea = sum(simplexvolume(node,tt(ring,:)));
            newq = simpqual(node,newt);
            if (abs(newarea-oldarea)/oldarea < 1e-3) && (newq>qi)
                tt(ring(1),:) = newt;
                q(ring(1)) = newq;
                tt(ring(2:end),:) = 0;
                elem(isoElem(ring(1)),:) = newt;
                elem(isoElem(ring(2:end)),1) = 0;
            else
                isIsoNode(isoNode(i)) = false;
            end
        else % linkpt contains 0, which means some triangles are deleted
            isIsoNode(isoNode(i)) = false;
        end
    end
end
%% Graph
if (nargin>2) && (gflag == 1)
    showmesh(node,elem);
    hold on; findnode(node,isIsoNode,'noindex','color','r')
end
%% Clean up
elem((elem(:,1) == 0),:) = [];
node(isIsoNode,:) = [];
indexMap = zeros(N,1);
indexMap(~isIsoNode)= 1:size(node,1);
elem = indexMap(elem);