function  vemPoisson3(cube, pde, option, varargin)

maxIt = option.maxIt;

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);
h = zeros(maxIt,1);
h(1) = option.h0;


% [node,elem, interfaceData] = interfacemesh3(cube, pde.phi, h(1));
% checkinterfacemesh3(node,elem,interfaceData)
% [face, face2elem] = getpolymesh(elem,interfaceData);
% clear elem interfaceData

%         [node, elem] = cubehexmesh(cube,h(k+1));
%         idx = (1:size(elem, 1))';
%         face =  [elem(:,[1,4,3,2]); elem(:,[6, 7, 8, 5]);...
%                  elem(:,[2,3,7,6]); elem(:,[1, 5, 8, 4]);...
%                  elem(:,[1,2,6,5]); elem(:,[4, 8, 7, 3])];
%         face2elem = [idx;idx;idx;idx;idx;idx];
%         clear elem;

[node, elem] = cubemesh(cube, h(1));
elem = fixorder3(node,elem);
idx = (1:size(elem,1))';
face = [elem(:,[2, 3, 4]);elem(:,[1, 4,3]);elem(:,[1,2,4]);elem(:,[1,3,2])];
face(:,4) = 0;
face2elem = [idx;idx;idx;idx];
clear elem

%% Vertual Element Method
for k = 1:maxIt
    N(k) = size(node,1);
    [u, info] = Poisson3VEM(node, face, face2elem, pde, option);
    errMax(k) = max(abs(u - pde.exactu(node)));
    if k < maxIt
          h(k+1) = h(k)/2;
          
%         [node, elem, interfaceData] = interfacemesh3(cube, pde.phi, h(k+1));
%         h(k+1)
%         checkinterfacemesh3(node,elem,interfaceData)
%         [face, face2elem] = getpolymesh(elem,interfaceData);
%         clear elem interfaceData

%         [node, elem] = cubehexmesh(cube,h(k+1));
%         idx = (1:size(elem, 1))';
%         face =  [elem(:,[1,4,3,2]); elem(:,[6, 7, 8, 5]);...
%                  elem(:,[2,3,7,6]); elem(:,[1, 5, 8, 4]);...
%                  elem(:,[1,2,6,5]); elem(:,[4, 8, 7, 3])];
%         face2elem = [idx;idx;idx;idx;idx;idx];
%         clear elem;
        
        [node, elem] = cubemesh(cube, h(k+1));
        elem = fixorder3(node,elem);
        idx = (1:size(elem,1))';
        face = [elem(:,[2, 3, 4]);elem(:,[1, 4, 3]);...
                elem(:,[1, 2, 4]);elem(:,[1, 3, 2])];
        face(:,4) = 0;
        face2elem = [idx;idx;idx;idx];
        clear elem
    end
end

set(gcf, 'Units', 'normal');
set(gcf, 'Position', [0.25, 0.25, 0.55, 0.4]);
showrate(1./h(1:k), errMax(1:k), 1, 'k-+', '||u-u_h||');


end
