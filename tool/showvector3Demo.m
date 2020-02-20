%% Examples of showvector3
clear;
help showvector3;

%% vector field
radius = @(p) sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
ufct = @(p)[3*p(:,1).*p(:,3)./radius(p).^2,...
            3*p(:,2).*p(:,3)./radius(p).^2,...
            (3*p(:,3).^2-radius(p).^2)./radius(p).^2];

[X, Y, Z] = ndgrid(-2:0.2:2);

% p is an N*3 array
p = [X(:) Y(:) Z(:)];
u = ufct(p);

%% Example 1: default input
figure(1);%#ok<*UNRCH>
h1 = showvector3(p,u);

%% Example 2: use quiver3 using the option structure
option.scale = 2.6;
option.plot = 'quiver3';
figure(2);
h2 = showvector3(p,u,option);
axis tight;

%% Example 3: modify parameters and figure properties
option.scale = 1;
option.plot = 'coneplot';
option.interp = 'linear';
option.stepSize = 0.4;
option.layer = 8;

figure(3);
h3 = showvector3(p,u,option,...
    'FaceAlpha',0.8,...
    'FaceColor',[0.5,0.8,0.8],...
    'EdgeAlpha',0.0);
h3.FaceLighting = 'gouraud';

% see Patch Properties and QuiverGroup Properties for more information
% http://www.mathworks.com/help/matlab/ref/patch_props.html
% http://www.mathworks.com/help/matlab/ref/quivergroupproperties.html

%% 
fullscreen = get(0,'ScreenSize');
set(gcf,'Position',[0 50 fullscreen(3)*0.5 fullscreen(4)*0.5])
set(gcf,'color','w')