function findnode(node,range,varargin)
%% FINDNODE highlights nodes in certain range.
%
%    FINDNODE(node,range) finds all nodes whose indices are in the range
%    array. It plots these nodes using larger black dots and display the indices.
%
%    FINDNODE(node) finds all nodes.
%   
%    FINDNODE(node,range,'noindex') skip the display of indices.
%
%    FINDNODE(node,range,'param','value','param','value'...) allows
%    additional patch param/value pairs to highlight the nodes.
%    
% Example
%     node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0];
%     elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];
%     subplot(1,2,1);
%     showmesh(node,elem);
%     findnode(node,[1 3],'noindex','color','r','MarkerSize',24)
%     subplot(1,2,2);
%     showmesh(node,elem);
%     findnode(node);
%
%   See also findelem, findnode3, findedge.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

hold on
dotColor = 'r.';
if (nargin==1) || (ischar(range) && strcmp(range,'all'))
    range = (1:size(node,1))'; 
    dotColor = 'k.';
end
if size(node,2) == 3 % 3-D case
    findnode3(node,range,varargin{:});
    return
end
% if isempty(range), return; end
if islogical(range) 
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end 
% if (nargin>2) && (strcmp(varargin{1},'coordinate'))
%     distance = (node(:,1)-range(1)).^2 + (node(:,2)-range(2)).^2;
%     [mindist,range] = min(distance); %#ok<ASGLU>
% end
h = plot(node(range,1),node(range,2),dotColor, 'MarkerSize', 20);
if nargin>2
    startidx = 1;
    if strcmp(varargin{1},'noindex') || strcmp(varargin{1},'index') || ...
       strcmp(varargin{1},'coordinate') || isnumeric(varargin{1})
        startidx = 2;
    end
    if nargin>3 && (strcmp(varargin{2},'noindex') || strcmp(varargin{2},'index'))
        startidx = 3;
    end
    if size(varargin,2)>=startidx
        set(h,varargin{startidx:end});
    end    
end
idxFlag = 1; % default show index number
if nargin>2 && isnumeric(varargin{1})
    shift = varargin{1};
else
    shift = [0.015 0.015];
end
if (nargin>3 && (strcmp(varargin{1},'noindex') || strcmp(varargin{2},'noindex'))) || ...
   (nargin>2 && (strcmp(varargin{1},'noindex')))
    idxFlag = 0;
end
if (nargin>3 && (strcmp(varargin{1},'index') || strcmp(varargin{2},'index'))) || ...
   (nargin>2 && (strcmp(varargin{1},'index')))        
    idxFlag = 1;
end
if idxFlag
   text(node(range,1)+shift(1),node(range,2)+shift(2),int2str(range), ...
     'FontSize',14,'FontWeight','bold');
end
hold off