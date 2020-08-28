%% An example of using bisect3 on unstructured 3D mesh in (0,1)^3
% The data is generated from delaunayTriangulation after R2013b
% The attribute names may be different before R2016a.
clear;close all;
N = randi([10,20]); 

%% randomly generating N points
box = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
coords = linspace(0.3,0.7,N);
n = floor(N/2);
node = [coords(randi(n,[1,n]))',coords(randi(n,[1,n]))',coords(randi(n,[1,n]))'];
node = [box; node];

%% Delaunay triangulation
dt = delaunayTriangulation(node);

%%
elem = dt.ConnectivityList;

%% initial grid
N0 = size(node,1);
HB = zeros(N0,4);
HB(1:N0,1:3) = repmat((1:N0)',1,3);

%% visualize

subplot(1,2,1);
h = tetramesh(elem,node,ones(size(elem,1),1));
set(h,'FaceAlpha',0.2,'FaceColor',[0.5,0.8,0.8]);    
view(3);
camproj('perspective');
axis off; axis equal; axis tight

%% refine using bisect3
NT = size(elem,1);
markedElem = randperm(NT, floor(NT/3)); % refine 1/3 of the elements
[node,elem,~,HB] = bisect3(node,elem,markedElem,[],HB);

%% visualize the finer mesh
subplot(1,2,2);
h = tetramesh(elem,node,ones(size(elem,1),1));
set(h,'FaceAlpha',0.1,'FaceColor',[0.5,0.8,0.8]);    
view(3);
camproj('perspective');
axis off; axis equal; axis tight
set(gcf,'color','w','Position', [100 200 800 300])
drawnow;