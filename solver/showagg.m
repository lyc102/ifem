function showagg(node,node2agg,agg2node,A,varargin)
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

N = size(node,1);
hold on
%
idx = find(node2agg~=0);
aggidx = (agg2node~=0);
% a simple prolongation
Nc = max(node2agg);
Pro = sparse(idx,node2agg(idx),1,N,Nc);
Ac = Pro'*A*Pro;  % graph of aggs
% coloring of aggregates
nodeC = node(agg2node(aggidx),:);
nodeSet = coloring(Ac);  % assign different color
% plot agg
showcoloring(Ac,nodeC,nodeSet);
colorstring = ['r' 'y' 'b' 'g' 'k' 'c' 'm' 'w'];
for k = 1:size(nodeSet,1)
    aggNumber = nodeSet{k};
    if ~isempty(aggNumber)
        [i,j] = find(Pro(:,aggNumber)); %#ok<*NASGU>
        findnode(node,i,'noindex','Color',colorstring(k),'MarkerSize',48);
    end
end
if nargin>4 && isnumeric(varargin{1})
    shift = varargin{1};
else
    shift = [0.015 0.015];
end
% text(node(idx,1)+shift(1),node(idx,2)+shift(2),int2str(node2agg(idx)), ...
%  'FontSize',14,'FontWeight','bold');
hold off
