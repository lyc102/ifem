function showcoloring(G,node,nodeSet)
%% SHOWCOLORING plot coloring vertices set
%
% Example
%
%     load lakemesh
%     N = size(node,1); NT = size(elem,1);
%     t2p = sparse([1:NT,1:NT,1:NT], elem(1:NT,:), 1, NT, N);
%     G = t2p'*t2p;
%     nodeSet = coloring(G);     
%     showcoloring(G,node,nodeSet);     
%
% See also: coloring, showgraphmatrix
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

showgraphmatrix(G,node);
colorstring = ['r' 'y' 'b' 'g' 'k' 'c' 'm' 'w'];
for k = 1:size(nodeSet,1)
    if any(nodeSet{k})
        findnode(node,nodeSet{k},'noindex','Color',colorstring(k));
    end
end