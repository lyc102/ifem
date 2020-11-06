function [nodeInterface, elemInterface, nodeBg, elemBg] = get2dpolymesh(node, elem, phi)
%% GET2DPOLYMESH gets the interface fitted polygonal mesh for a 2d quadtree mesh


%% Prepare
if ~iscell(elem); elem = num2cell(elem,2); end

elemBg = elem;
nodeBg = node;
T = auxstructurepoly(elem);
edge = T.edge;
edge2elem = double(T.edge2elem);

N = size(node,1);
NT = size(elem,1);

%% compute the level set value at each vertex
% 
phiValue = phi(node);
vSign = msign(phiValue);

%% Step 1: Find points ON interface 
% Find the intersection points between edges and the interface
isCutEdge = (vSign(edge(:,1)).*vSign(edge(:,2))<0);
A = node(edge(isCutEdge,1),:);
B = node(edge(isCutEdge,2),:);
nodeCut = findintersectbisect(phi,A,B);
Ncut = size(nodeCut, 1);
ixCutNode = N+(1:Ncut)';
ixIntersectNode = find(abs(vSign) <= eps);
Nintersect = size(ixIntersectNode,1);
vSign(N+1:N+Ncut) = 0;
vSign(N+Ncut+1:N+Ncut+Nintersect) = 0;

%% Step 2: find interface elem and nodes
isInterfaceElem = false(NT,1);  
interfaceElem = edge2elem(isCutEdge,[1,2]);
isInterfaceElem(interfaceElem) = true;
% new vertices needed to be merged into interface element
elemUpdate = accumarray(interfaceElem(:),[ixCutNode;ixCutNode], [NT,1],@(x){x'});
elem = cellfun(@horzcat, elem, elemUpdate, 'UniformOutput', false);
elemVertNum = cellfun('length',elem);

vSign2elem = mat2cell(vSign([elem{:}]),elemVertNum);
vSign2elem = cellfun(@transpose, vSign2elem, 'UniformOutput', 0);
isInterfaceElem(cellfun(@(x) sum(x==0)==2, vSign2elem)) = true; % 2 vertices on interface
ixIFElem = find(isInterfaceElem);
% elemOnInterface = elem(isInterfaceElem,:);

nodeInterface = [node; nodeCut];

%% cut polygonal element into sub elements for computation
numInterfaceElem = length(ixIFElem);
elem(NT+1:NT+numInterfaceElem) = cell(numInterfaceElem,1);


for i = 1:numInterfaceElem
   ix = ixIFElem(i);
   ixLocElem1= [find(vSign2elem{ix}==0), find(vSign2elem{ix}==1)];
   ixLocElem2 = [find(vSign2elem{ix}==0), find(vSign2elem{ix}==-1)];
%    elem{ix} = {elem{ix}(ixLocElem1); elem{ix}(ixLocElem2)};
    elem{NT+i} = elem{ix}(ixLocElem2);
    elem{ix} = elem{ix}(ixLocElem1);
end

% ixIFElemNew = [ixIFElem; ((NT+1):(NT+numInterfaceElem))'];
elemInterface = fixorientationpoly(nodeInterface,elem);

end