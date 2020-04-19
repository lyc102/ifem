function movingcircle
%% MOVINGCIRCLE will track the interface defined by x^2+y^2 = (0.5-t)^2
%--------------------------------------------------------------------------
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details. 
%--------------------------------------------------------------------------

%%
close all; clear variables;
%---------------------- Parameters ----------------------------------------
figure(1); set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.8,0.4]);
t = 0; dt = 0.025; maxIt = 120;
%---------------------- Initial Grid --------------------------------------
node = [-1 -1; 1 -1; 1 1; -1 1];
elem = [2 3 1; 4 1 3];
% N0 = size(elem,1);
for i = 1:6
	[node,elem] = uniformbisect(node,elem);
end
u = -sign(f(node,t));
subplot(1,2,1); showmesh(node,elem); pause(0.1)
subplot(1,2,2); showsolution(node,elem,u,[0,90]); 
for k = 1:maxIt
    % move interface every four iterations	
    if mod(k,4) == 0
        t = t + dt; 
    end
	%---------- detect element cross interface or away from interface -----
	eta = abs(sign(f(node(elem(:,1),:),t)) + sign(f(node(elem(:,2),:),t))...
            + sign(f(node(elem(:,3),:),t)));
	refineElem = find(eta < 3);
    coarsenElem = find(eta == 3);
    %---------- refine elements cross the interface -----------------------
    [node,elem] = bisect(node,elem,refineElem);
    u = -sign(f(node,t));
    subplot(1,2,1);  showmesh(node,elem); pause(0.0125)
    subplot(1,2,2);  showsolution(node,elem,u); view(2);
    %---------- coarsen elements away from the interface ------------------
    [node,elem] = coarsen(node,elem,coarsenElem);
    u = -sign(f(node,t));
    %---------- plot mesh and function ------------------------------------
    subplot(1,2,1);  showmesh(node,elem); pause(0.0125)
    subplot(1,2,2);  showsolution(node,elem,u); view(2);
end
end %%  ------- End of MOVINGCIRCLE --------------------------------

%% -------------------- Sub functions called by MOVINGCIRCLE ---------------
function z = f(p,t)
z = sum(p.^2,2) - (0.85-t)^2;
end