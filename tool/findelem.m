function findelem(node,elem,range,varargin)
%% FINDELEM highlights some elements
%
%    FINDELEM(node,elem,range) finds all elements whose indices are in the
%    range array by displaying these elements in yellow. The size of the marker
%    is proportional to the number of digits of the indices. The font will
%    magnify as one resizes the figure.
%
%    FINDELEM(node,elem) finds all elements.
%   
%    FINDELEM(node,elem,range,'noindex') skip the display of indices.
%
%    FINDELEM(node,elem,range,'param','value','param','value'...) allows
%    additional patch param/value pairs to highlight the elements.
%    
% Example:
%     node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0];
%     elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];
%     subplot(1,2,1);
%     showmesh(node,elem);
%     findelem(node,elem,1,'index','FaceColor','r');
%     subplot(1,2,2);
%     showmesh(node,elem);
%     findelem(node,elem);
%
%   See also findelem3, findnode3, findedge.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on

if (nargin==2) || isempty(range) || (ischar(range) && strcmp(range,'all'))
    range = (1:size(elem,1))'; 
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
nV = size(elem,2);
if nV == 3
    center=(node(elem(range,1),:)+node(elem(range,2),:)+node(elem(range,3),:))/3;
elseif nV == 4
    center=(node(elem(range,1),:)+node(elem(range,2),:)+node(elem(range,3),:)+node(elem(range,4),:))/4;
end
if length(range) < size(elem,1)
    x = reshape(node(elem(range,:),1),size(range,1), size(elem,2))';
    y = reshape(node(elem(range,:),2),size(range,1), size(elem,2))';
    h = patch(x,y,'y');
    if nargin > 3
        if strcmp(varargin{1},'noindex') || strcmp(varargin{1},'index')
           if size(varargin,2)>=2
                set(h,varargin{2:end});
           end
        else
            set(h,varargin{:});
        end
    end
end
nDigit =ceil(log10(range+1));
oSize = 400*nDigit;
f = gcf;
set(f, 'Units', 'pixel');
fHeight = f.Position(4);
if (nargin <=3) || ~(strcmp(varargin{1},'noindex'))
    if size(node,2) == 2
        s = scatter(center(:,1),center(:,2),oSize,'o','LineWidth',1,'MarkerEdgeColor','k',...
            'MarkerFaceColor','y');
        if max(f.Position(3:4)) <= 1
            fSz = 20*f.Position(4); % font size is proportional to height
        else
            fSz = log(f.Position(3)*f.Position(4));
        end
        tShift = 0.02*(max(nDigit) + 1 - nDigit);
        tShift(nDigit==1) = 0.5*tShift(nDigit==1);
        t = text(center(:,1)-tShift,center(:,2),int2str(range),...
            'Fontsize',fSz,'FontWeight','bold','Color','k');
        nElemShown = length(t);
        set(f,'SizeChangedFcn',@resizeCallback);
        set(f,'Visible','on');
        drawnow;
    elseif size(node,2) == 3 % surface mesh
        scatter3(center(:,1),center(:,2),center(:,3),oSize,'o','LineWidth',1,...
            'MarkerEdgeColor','k','MarkerFaceColor','y');    
        text(center(:,1)-0.02*nDigit,center(:,2),center(:,3),int2str(range),'FontSize',12,...
        'FontWeight','bold','Color','k');        
    end
end
hold off
%% relative position of numbers
% add by Shuhao Cao
    function resizeCallback(f, ~)
        if max(f.Position(3:4)) <= 1 % relative
            newFSz = 80*f.Position(4); % font size is proportional to height
            newOSize = 10*oSize*(f.Position(4));
            newTshift = 4*tShift*f.Position(4);
        else % pixel
            newFSz = log(f.Position(3)*f.Position(4));
            newOSize = oSize*(f.Position(4)/fHeight);
            newTshift = tShift/(f.Position(4)/fHeight);
        end
        % change font size accordingly
        set(t,'FontSize',newFSz);
        set(s,'SizeData',newOSize);
        for i = 1:nElemShown
            set(t(i), 'Position', [center(i,1)-newTshift(i), center(i,2), 0]);
        end
    end

end