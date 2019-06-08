function nodeSet = coloring(G,nc,qp)
%% COLORING color vertices of a graph
%
% nodeSet = coloring(G,nc) color the vertices of the graph G by nc colors.
% Vertices in the same color are disjointed except the last one. The input
% graph G is represented by the sparse pattern of a symmetric sparse
% matrix. The values of G is not used. The output nodeSet is a cell
% array and nodeSet{k} is the indices of vertices in the k-th color. By
% default, nc = 6 if not specified by the input.
%
% nodeSet = coloring(G,nc,qp) chose the vertices with a larger value qp
% first. For example, qp can be the degree of vertices. Then nodeSet{1}
% will pick up nodes with a local max degree. If no qp, random values are
% used.
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
% See also: showcoloring, showgraphmatrix 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Set up
N = size(G,1);
N0 = min(floor(N/100),20);  % number of the coarest nodes
if ~exist('nc','var'), nc = 8; end
nodeSet = cell(nc,1);
if ~exist('qp','var')
    qp = rand(N,1); 
else
    idx = (qp ~= 0);
    qp(idx) = (1 + 0.1*rand(sum(idx),1)).*qp(idx);
end

%% Coloring
for k = 1:nc-1
    isC = false(N,1);       % C: coarse nodes
    isF = false(N,1);       % F: fine nodes
    isU = (qp>0);           % U: undecided nodes
    if ~any(isU)            % no more vertices left
        break
    end
    for m = 1:round(log10(N))
        % Mark all available nodes
        isS = isU;          % selected all undecided vertices
        S = find(isS); 

        % Find marked nodes with local maximum quality
        [i,j] = find(triu(G(S,S),1));   % i,j and i<j: edges of subgraph S
        idx = (qp(S(i)) >= qp(S(j)));   % compare quality of connected vertices
        isS(S(j(idx))) = false;         % remove vertices with smaller qp
        isS(S(i(~idx))) = false;        % if qualities are equal, keep the nodes with smaller index
        isC(isS) = true;                % set selected nodes as coarse nodes

        % Remove coarse nodes and neighboring nodes from undecided set
        [i,j] = find(G(:,isC)); %#ok<*NASGU>
        isF(i) = true;                  % neighbor of C nodes are F nodes
        isU(isF | isC) = false;         % remove current C and F from U
        if sum(isU) <= N0               
            break;    % only break the inner loop
        end
    end
    nodeSet{k} = find(isC);         
    qp(nodeSet{k}) = 0;    
end
nodeSet{k+1} = find(qp>0);  % all left vertices are added to the last color
% nodeSet = nodeSet{1:k+1};