function showgraphmatrix(G,node)
%% SHOWGRAPHMATRIX displays a planar graph
%
%    showgraphmatrix(G,node) displays a planar undirected graph represented
%    by a sparse matrix G.
%
%   See also showgraph, showmesh, findedge
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[i,j] = find(triu(G,1));
showgraph(node,[i j]);