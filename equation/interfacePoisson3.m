function [u, w, AE, AI, info] = interfacePoisson3(node,elem, pde,interfaceData,option)
%% INTERFACEPOISSON3 Poisson equation: P1 linear element in Virtual Element Methods.
%
% Output: AE is the assembling matrix for exterior elements
%         AI is the assembling matrix for interior elements
%   u = INTERFACEPOISSON3(node,elem,pde,E) produces the linear finite element
%   approximation of the interface Poisson equation
% 
%       -div(\beta*grad(u))=f  in \Omega\ \Gamma, with 
%       Dirichlet boundary condition:  u=g_D on \partial\Omega, 
%       Jump conditions across the interface:  [u]_{\Gamma}=u^+-u^- = q_0 on \Gamma,
%       [\beta u_n]_{\Gamma}= \beta^+ u_n^+ -\beta^+- u_n^- =q_1 on \Gamma
%
% See also: interfacePoisson

elemIdx = interfaceData.elemIdx;
face = interfaceData.face;
face2elem = interfaceData.face2elem;
vSign = interfaceData.vSign;
interface = interfaceData.interface;

if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
%% classfication for elements
N = size(node,1);
NC = size(elem,1);
% uncut elements
oldElem = elem(~(elemIdx==0),:);
idx = find(elemIdx~=0);

face = [face;...
        oldElem(:,1) oldElem(:,4) oldElem(:,3) oldElem(:,2);...
        oldElem(:,6) oldElem(:,7) oldElem(:,8) oldElem(:,5);...
        oldElem(:,2) oldElem(:,3) oldElem(:,7) oldElem(:,6);...
        oldElem(:,1) oldElem(:,5) oldElem(:,8) oldElem(:,4);...
        oldElem(:,1) oldElem(:,2) oldElem(:,6) oldElem(:,5);...
        oldElem(:,4) oldElem(:,8) oldElem(:,7) oldElem(:,3)];
%NF = size(face,1);
face2elem = [face2elem;idx;idx;idx;idx;idx;idx];

isTriFace = (face(:, 4) == 0);
MN = max(face2elem);

elem2node = sparse(face2elem(isTriFace)*ones(1,3), face(isTriFace,1:3), 1, MN, N)...
    + sparse(face2elem(~isTriFace)*ones(1,4), face(~isTriFace,1:4), 1, MN, N);
%elem2face = sparse(face2elem, (1:NF)', 1, MN, NF);

elem2node = (elem2node > 0);
NV = elem2node*ones(N,1);
centroid = elem2node*node./[NV, NV, NV];
NVS = unique(NV);

coefficient = 1/3*ones(size(face,1),1);
coefficient(isTriFace) = 1/6;
normal = faceNormal(node,face); % face normal vectors 
volume = accumarray(face2elem, coefficient.*dot(node(face(:,1),1:3),normal(:,1:3),2)); % why?
h = volume.^(1/3);
B = 1/6*normal./(h.^2*ones(1, 3));
B(~isTriFace,:) = 3/2 * B(~isTriFace,:);

for nv = NVS
    isCurrentFace = (NV(face2elem) == nv);
    isCurrentElem = false(MN, 1);
    isCurrentElem(face2elem(isCurrentFace)) = true;
    currentElem2node = elem2node(isCurrentElem,:);
    [N1, N2] = size(currentElem2node,1);
    [~, J] = find(currentElem2node');
    currentElem = reshape(J, nv, [])';
    currentElem2localIdx = sparse(repmat((1:N1)',1, Nv),currentElem,ones(N1,1)*(1:nv), N1, N2);
    
    isCTFace = isTriFace & isCurrentFace;
    N1 = sum(isCTFace);
    tFace = face(isCTFace,1:3);
    subs1 = repmat([ones(N1, 1); 2*ones(N1,1); 3*ones(N1,1)], 1,3);   
    subs3 = repmat(face2elem(isCTFace),1, 9);
    subs2 = [repmat(tFace(:,1),1, 3);repmat(tFace(:,2),1, 3);repmat(tFace(:,3),1, 3)];
    subs2 = currentElem2localIdx((subs2-1)*N1+ subs3);
    val = repmat([B(isCTFace,1);B(isCTFace,2);B(isCTFace,3)], 1, 3);
    
    BP = accumarray({subs1, subs2, subs3}, val);
    
    isCSFace = ~isTriFace & isCurrentFace;
    N1 = sum(isCSFace);
    sFace = face(isCSFace,:);
    subs1 = repmat([ones(N1, 1); 2*ones(N1,1); 3*ones(N1,1)], 1,3);
    subs3 = repmat(face2elem(isCTFace),1, 12);
    subs2 = [repmat(sFace(:,1),1, 3);repmat(sFace(:,2),1, 3);repmat(sFace(:,3),1, 3);repmat(sFace(:,4),1, 3)];
    subs2 = currentElem2localIdx((subs2 - 1)*N1 + subs3);
    val = repmat([B(isCSFace,1);B(isCSFace,2);B(isCSFace,3)], 1, 4);
    BP = BP + accumarray({subs1, subs2, subs3}, val);      
end


end
