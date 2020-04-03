function [h1, h2] = showboundary3(node,elem,expr,varargin)
%% SHOWBOUNDARY3 boundary surface mesh plot
%
%    showboundary3(node,elem) displays the boundary surface mesh of a
%    3-dimensional tetrahderal mesh given by node and elem matrices. The
%    boundary is found by the function findboundary3.
%
%    showboundary3(node,elem,expr) displays the boundary of parts of the
%    mesh specificed by the expression. For example,
%    showboundary3(node,elem,'tempvar(x>=0 & y>=0 & z>=0)') only shows the boundary
%    mesh for tetrahedra not in the first quadrant. 
%
%    showboundary3(node,elem,expr,viewangle) changes the display angle. The
%    default view angle is view(3).
%
%    showboundary3(node,elem,expr,'param','value','param','value'...)
%    allows additional patch param/value pairs to be used when displaying
%    the mesh. For example, the default transparency parameter is set to
%    0.5. You can overwrite this value by using the param pair
%    ('FaceAlpha', value). The value has to be a number between 0 and 1. 
%    For all other possible parameters, see Patch Properties.
%    https://www.mathworks.com/help/matlab/ref/matlab.graphics.primitive.patch-properties.html
%    Addition param: 'Output', when output is set to 1, a vector formatted
%    pdf is printed to the current Matlab directory.
%   
%    To display the full tetrahedral mesh in 3-d, use showmesh3. Notice that the
%    3-d graphical visualization is slow for large mesh data. 
%
%   Example:
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     [node,elem] = uniformbisect3(node,elem);
%     subplot(1,3,1);
%     showboundary3(node,elem);
%     subplot(1,3,2);
%     showboundary3(node,elem,'~(x>=0 & y>=0 & z>=0)',[59,20]); 
%     subplot(1,3,3);
%     showboundary3(node,elem,'~(x>=0 & y>=0 & z>=0)',[59,20], 'Output', 1); 
%
%   See also showboundary3, showsolution3, showmesh3.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.
%
% Modified from simpplot.m in distmesh by Per-Olof Persson.

%% Find boundary faces and cut faces
[tempvar,allBdFace] = findboundary3(elem); %#ok<*ASGLU>
if (nargin >= 3) && (any(expr))
    x = node(:,1);  y = node(:,2);  z = node(:,3); %#ok<*NASGU>
    incl = find(eval(expr));
    elem = elem(any(ismember(elem,incl),2),:);
end
[tempvar, bdFace] = findboundary3(elem);
cutFace = setdiff(bdFace,allBdFace,'rows');
%% Plot surface mesh
% bcol=.9*ones(1,3); icol=.6*ones(1,3);
bcol=[0.45 1 0.45]; icol=.75*bcol;
h1 = showmesh(node,bdFace,'facecolor',bcol,'edgecolor','k');
hold on
if ~isempty(cutFace)
    h2 = showmesh(node,cutFace,'facecolor',icol,'edgecolor','k');
end
isOutput = false;
if ~isempty(varargin)
    outputOption = cellfun(@isequal,varargin,repmat({'output'},size(varargin)));
    
    if any(outputOption)
        if varargin{find(outputOption)+1}
            isOutput = true;
            idxisOutput = find(outputOption);
            varargin([idxisOutput, idxisOutput+1]) = [];
        else
            idxisOutput = find(outputOption);
            varargin([idxisOutput, idxisOutput+1]) = [];
        end
    end
    
    if isnumeric(varargin{1})
        view(varargin{1});
        if nargin>5
            set(h1,varargin{2:end});
            if ~isempty(cutFace)
                set(h2,varargin{2:end});
            end
        end
    else
        set(h1,varargin{1:end});
        if ~isempty(cutFace)
            set(h2,varargin{1:end});
        end
    end
end
% light;lighting phong;
cameratoolbar;
drawnow;

if isOutput; printOutput; end

function printOutput
    %% shrink the axes to figure
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    %% output size
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    
    %% print to pdf
    filename = strcat('tetramesh_',...
        '[', num2str(min(node(:,1))), ',', num2str(max(node(:,1))), ']',...
        'x','[', num2str(min(node(:,2))),',', num2str(max(node(:,2))),']',...
        'x', '[', num2str(min(node(:,3))),',', num2str(max(node(:,3))),']',...
        '_node',num2str(size(node,1)),...
        '_elem',num2str(size(elem,1)));
    print(fig,filename,'-dpdf','-painters');
end

end