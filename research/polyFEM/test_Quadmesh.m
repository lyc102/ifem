%% Testing ground for quad mesh adaptive refinement 
% This is the testing utilizing VEM's polytopal structure
% without using the tree-structure to maintain the mesh data.
%
% Outer loop: # of vertices
% inner loop:
%
% Step 1: find the center/centroid P of marked elements
% Remark: current version using naive way of locating centroid
%         
% Step 2: quad-sect using P, find duplicate nodes
% Step 3: form new elements replacing the marked elements using cell
% Step 4: add vertices to the neighbor elements through edges
%
% Two sets of mesh are kept:
% 1. For refinement elem: the first 4 vertices are always the square clockwise,
% after those are the extra vertices for polygon.
% 2. For computation elemRefine: re-ordered vertices mesh.

clear; close all;

%% 
[node,elemOrig] = squarequadmesh([0,1,0,1],0.25);
figure(1); 
showmesh(node,elemOrig);
findnode(node);
set(gcf,'color','w','Position', [100 600 500 500])

% local ordering should always be kept as
% 2 -- 3
% |    |
% 1 -- 4

T = auxstructurequad(elemOrig);

elem = num2cell(elemOrig,2); % convert elem matrix to cell

elemVertexNumber = cellfun('length',elem);

N = size(node,1);


%% 
% Nv = min(elemVertexNumber):max(elemVertexNumber)
Nv=4;
    
% find polygons with Nv vertices
idx = find(elemVertexNumber == Nv); % index of elements with Nv vertices
NT = length(idx); % the number of elements having the same number of vertices
% vertex index and coordinates
vertex = cell2mat(elem(idx));

x = reshape(node(vertex,1),NT,Nv);
y = reshape(node(vertex,2),NT,Nv);

xmax = max(x,[],2); % the rightmost vertices' x-coords for each element
ymax = max(y,[],2); % the uppermost vertices' y-coords for each element
xmin = min(x,[],2); % the leftmost vertices' x-coords for each element
ymin = min(y,[],2); % the lowermost vertices' y-coords for each element

xc = 0.5*(xmax + xmin);
yc = 0.5*(ymax + ymin);

x_nb = circshift(x,[0,-1]);
y_nb = circshift(y,[0,-1]);

%% marking elements
idxMarkedElem = idx(randi(NT, [floor(NT/3), 1]));
% idxMarkedElem = idx([1,3]);
NmarkedElem = length(idxMarkedElem);

NnodeNew = 5*NmarkedElem; 
% number of the new nodes per Nv consisting duplicated nodes


% o -- 2 -- o
% |    |    |
% 1 -- 5 -- 3
% |    |    |
% o -- 4 -- o
% o: original vertices
% 1 to 5: local indices of newly added vertices
% the order is consistent with neighbor ordering

nodeNew(1:5:NnodeNew,:) = [xmin(idxMarkedElem,1), yc(idxMarkedElem,1)]; % left
nodeNew(2:5:NnodeNew,:) = [xc(idxMarkedElem,1), ymax(idxMarkedElem,1)]; % top
nodeNew(3:5:NnodeNew,:) = [xmax(idxMarkedElem,1), yc(idxMarkedElem,1)]; % right
nodeNew(4:5:NnodeNew,:) = [xc(idxMarkedElem,1), ymin(idxMarkedElem,1)]; % bottom
nodeNew(5:5:NnodeNew,:) = [xc(idxMarkedElem,1), yc(idxMarkedElem,1)]; % center

[nodeNewU, idxN, idxU] = unique(nodeNew,'rows','legacy');
% nodeNewU = nodeNew(idxN,:) and nodeNew = nodeNewU(idxU,:)

idxNodeNew = N+1:N+NnodeNew; % unique indices of new nodes


%% generate new elements on marked elements
elem2nodeNew = reshape(idxNodeNew(idxU), [5, NmarkedElem])';
% elem2nodeNew(i,:): i-th marked elem's newly added vertices
% elem2nodeNew(:,1): center of the marked elements
% elem2nodeNew(:,2:5): counterclockwisely ordered new edge vertices from
% the left edge locally

node = [node; nodeNewU];

elemMarked = cell2mat(elem(idxMarkedElem));

% localNewVertex = [2 1 3; 5 1 2; 4 1 5; 3 1 4]; %not very helpful

elemNew1 = horzcat(elemMarked(:,1), elem2nodeNew(:,[1 5 4]));
elemNew2 = horzcat(elem2nodeNew(:,1), elemMarked(:,2), elem2nodeNew(:,[2 5]));
elemNew3 = horzcat(elem2nodeNew(:,[5 2]), elemMarked(:,3), elem2nodeNew(:,3));
elemNew4 = horzcat(elem2nodeNew(:,[4 5 3]), elemMarked(:,4));

NelemNew = 4*NmarkedElem;

%% change elements' vertices on unrefined neighbor elements
% no need to check whether a vertex is added twice b/c
% the marked elements will get over-written in the next routine

elemMarkedNB = T.neighbor(idxMarkedElem,:);
% idxNb = ~(elemMarkedNB==idxMarkedElem);
elemMarkedNB = num2cell(elemMarkedNB,2);

for i = 1:size(elemMarked,1)
    for j = 1:size(elemMarkedNB{i},2)
        elem{elemMarkedNB{i}(:,j)} = ...
         [elem{elemMarkedNB{i}(:,j)}, elem2nodeNew(i,j)];
    end
end

%% vertically concat all new elements into a cell
elemNew = reshape([elemNew1'; elemNew2'; elemNew3'; elemNew4'], [4, NelemNew])';

elem(idxMarkedElem) = mat2cell(elemNew, repmat(4,[1,NmarkedElem]));

% elemTemp = vertcat(elem{:}); 
% above doesn't work if elem is of different size when nested

outTemp = cellfun(@(x) num2cell(cat(2,x),2),elem,'UniformOutput',false);
elemRefine = cellflat(outTemp)';

% celldisp(elemRefine)

%% reorder the polygon with vertices > 4
% to-do: rewrite the for loop to cellfun
elemVertexNumberRefine = cellfun('length',elemRefine);
idxReorder = find(elemVertexNumberRefine > 4); 

for i = 1:length(idxReorder)
    % compute the center (simple mean suffices) and sorting atan2
    elemTemp = elemRefine{idxReorder(i)};
    center = mean(node(elemTemp,:),1);
    vert2center = node(elemTemp,:) - center;
    [~, idxOrdered] = sort(atan2(vert2center(:,1), vert2center(:,2)));
    elemRefine{idxReorder(i)} = elemTemp(idxOrdered);
end


Nrefine = size(node,1); % update number of nodes

%%
figure(2);
showmeshpoly(node,elemRefine);
plot(xc,yc,'o');
findnode(node);
set(gcf,'color','w','Position', [650 600 500 500])
