%% MOVINGSPHERE
%
% The mesh will track the interface defined by x^2+y^2+z^2 = (0.75-t)^2
%
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details. 

close all; clear variables;
%% Parameters
figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.8,0.4]);
t = 0; dt = 0.1; maxIt = 16;

%% Generate an initial mesh 
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
for k = 1:2
    [node,elem,~,HB] = uniformbisect3(node,elem,[],HB);
end

%% Adaptive tracking of the moving interface
for k = 1:maxIt
    % move interface every four iterations	
    if mod(k,4) == 0
        t = t + dt; 
    end
	% detect element cross interface or away from interface
	eta = abs(sign(f(node(elem(:,1),:),t)) + sign(f(node(elem(:,2),:),t))...
            + sign(f(node(elem(:,3),:),t)) + sign(f(node(elem(:,4),:),t)));
	refineElem = find(eta < 4);
    coarsenElem = find(eta == 4);
    % refine elements cross the interface
    [node,elem,~,HB] = bisect3(node,elem,refineElem,[],HB);
    u = -sign(f(node,t));
    subplot(1,2,1); 
    showboundary3(node,elem,'~(x<=0 & y<=0)');
    pause(0.05)
    subplot(1,2,2);  
    showsolution3(node,elem,u,'~(x<=0 & y<=0)'); 
    colorbar;
    % coarsen elements away from the interface
    [node,elem,HB] = coarsen3(node,elem,coarsenElem,HB);
    u = -sign(f(node,t));
    subplot(1,2,1); 
    showboundary3(node,elem,'~(x<=0 & y<=0)');
    pause(0.05)
    subplot(1,2,2);  
    showsolution3(node,elem,u,'~(x<=0 & y<=0)'); 
    colorbar;
end
%% -------------------- Sub functions called by MOVINGCIRCLE ---------------
function s = f(p,t)
    s = sum(p.^2,2) - (0.75-t)^2;
end