%% MOVINGSPHERE
%
% The mesh will track the interface defined by x^2+y^2+z^2 = (0.75-t)^2
%
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details. 

close all; clear variables;
%% Parameters
figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.8,0.4]);
t = 0; dt = 0.025; maxIt = 50;

%% Generate an initial mesh 
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
[node,elem,HB] = delmesh(node,elem,'x<0 & y<0',HB);
bdFlag = setboundary3(node,elem,'Neumann');
[node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB);
% [node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB);

%% Adaptive tracking of the moving interface
for k = 0:maxIt
    % move interface every four iterations	
    if mod(k,2) == 0 && k > 4
        t = t + dt; 
    end
	% detect element cross interface or away from interface
	eta = abs(sign(f(node(elem(:,1),:),t)) + sign(f(node(elem(:,2),:),t))...
            + sign(f(node(elem(:,3),:),t)) + sign(f(node(elem(:,4),:),t)));
	refineElem = find(eta < 4);
    coarsenElem = find(eta == 4);
    % refine elements cross the interface
    [node,elem,bdFlag,HB] = bisect3(node,elem,refineElem,bdFlag,HB);
    u = -sign(f(node,t));
    subplot(1,2,1); 
    [tempvar,bdFace] = findboundary3(elem,bdFlag);
    showmesh(node,bdFace); 
    pause(0.0125)
    subplot(1,2,2);  
    showsolution(node,bdFace,u); 
    colorbar;
    % coarsen elements away from the interface
    [node,elem,bdFlag,HB] = coarsen3(node,elem,coarsenElem,bdFlag,HB);
    [tempvar,bdFace] = findboundary3(elem,bdFlag);
    u = -sign(f(node,t));
    subplot(1,2,1); 
    showmesh(node,bdFace); 
    pause(0.0125)
    subplot(1,2,2);  
    showsolution(node,bdFace,u); 
    colorbar;
end
%% -------------------- Sub functions called by MOVINGSPHERE ---------------
function s = f(p,t)
    s = sum(p.^2,2) - (0.75-t)^2;
end