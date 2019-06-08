function  vemInterfacePoisson3(cube, pde, option, varargin)

maxIt = option.maxIt;
h = zeros(maxIt,1);
h(1) = option.h0;

[node,elem, interfaceData] = interfacemesh3(cube, pde.phi, h(1));

%% Initialize err
errL2 = zeros(maxIt,1);   errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1); errMax = zeros(maxIt,1);
errTime = zeros(maxIt,1); solverTime = zeros(maxIt,1); 
assembleTime = zeros(maxIt,1); meshTime = zeros(maxIt,1); 
itStep = zeros(maxIt,1);  stopErr = zeros(maxIt,1); flag = zeros(maxIt,1);
N = zeros(maxIt,1);

%% Vertual Element Method
for k = 1:maxIt
    N(k) = size(node,1);
    [uEh, w, AE, AI, info] = interfacePoisson3VEM(node, elem, pde, interfaceData, option);
    interface = interfaceData.interface;
    elemIdx = interfaceData.elemIdx;

    isInterfaceNode = false(N(k),1);
    isInterfaceNode(interface(:)) = true;
    isExteriorNode = false(N(k), 1);
    isExteriorNode(elem(elemIdx == 1,:)) = true;
    isExteriorNode(isInterfaceNode) = true;
    isInteriorNode = false(N(k), 1);
    isInteriorNode(elem(elemIdx == -1,:)) = true;
    isInteriorNode(isInterfaceNode) = true;
    uIh = uEh;
    uIh(isInteriorNode) = uIh(isInteriorNode) - w(isInteriorNode);
    uEI = zeros(N(k), 1);
    uII = zeros(N(k), 1);
    uEI(isExteriorNode) = pde.exactuplus(node(isExteriorNode,:));
    uII(isInteriorNode) = pde.exactuminus(node(isInteriorNode,:));
    erruIuh(k) = sqrt((uEh - uEI)'*AE*(uEh - uEI) + (uIh - uII)'*AI*(uIh - uII));
    errMax(k) = max(max(abs(uEh(isExteriorNode) - uEI(isExteriorNode))),...
                    max(abs(uIh(isInteriorNode) - uII(isInteriorNode))));
    if k < maxIt
         h(k+1) = h(k)/2;
        [node, elem, interfaceData] = interfacemesh3(cube, pde.phi,h(k+1));
    end
end

set(gcf, 'Units', 'normal');
set(gcf, 'Position', [0.25, 0.25, 0.55, 0.4]);
showrate2(N(1:k), erruIuh(1:k), 1, 'b-*', '||u-u_h||_A',...
          N(1:k),errMax(1:k), 1, 'k-+', '||u-u_h||_{\infty}');
