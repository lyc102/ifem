function [u, info] = Poisson3VEM(node, face, face2elem, pde, option)
%% POISSONVEM3 Solve Poisson equation by linear vem.
%
%   u = PoissonVEM3(node, face, face2elem, pde) produces the linear vertual
%   element apprximation of the Poisson equation
%
%   The mesh is a polyhedra mesh which is given by face and face2elem. Here
%   we assume all of the faces are triangles or all of them are squares.


% number of nodes and povlyhedra
N = size(node, 1); NP = max(face2elem);
NF = size(face ,1);

Ndof = N;

isTriFace = face(:,4) == 0;
elem2node = sparse(face2elem(isTriFace)*ones(1,3),face(isTriFace,1:3),1, NP, N)...
 + sparse(face2elem(~isTriFace)*ones(1,4), face(~isTriFace,1:4),1, NP, N);

elem2node = (elem2node > 0);
NV = elem2node*ones(N,1);
centroid = elem2node*node./[NV, NV, NV];


coefficient = 1/6*ones(NF, 1);
coefficient(~isTriFace) = 2*coefficient(~isTriFace);
normal = facenormal(node,face); % face normal vectors 

%% \int_V div(x, y, z)^T = int_S (x,y,z)n 
volume = accumarray(face2elem, coefficient.*dot(node(face(:,1),:),normal,2)); % divergence therorem
h = volume.^(1/3);

B = 1/6*normal./(h(face2elem).^2*ones(1, 3));
B(~isTriFace,:) = 3/2*B(~isTriFace,:);

nnz = sum(NV.^2);
ii = zeros(nnz,1);
jj = zeros(nnz,1);
ss = zeros(nnz,1);
index = 0;
b = zeros(Ndof,1);

uNV = unique(NV);
for kk = 1:length(uNV)
    nv = uNV(kk);
    % Deal with elem with nv vertices
    isCurrentFace = (NV(face2elem) == nv);
    isCurrentElem = false(NP, 1);
    isCurrentElem(face2elem(isCurrentFace)) = true;
    
    % Current elem to node matrix
    currentElem2node = elem2node(isCurrentElem, :);
    CNP = size(currentElem2node, 1);
    currentElemLocalIdx = zeros(NP, 1);
    currentElemLocalIdx(isCurrentElem) = 1:CNP;
    BP = zeros(CNP,3,nv);

    [I, ~] = find(currentElem2node');
    currentElem = reshape(I, nv, [])';
    currentNodelocalIdx = sparse(repmat((1:CNP)', 1, nv), currentElem, ones(CNP, 1)*(1:nv), CNP, N);
    % Deal with current triangle face cases
    isCTFace = isCurrentFace & isTriFace;
    NN = sum(isCTFace);
    if NN > 0
        tFace = face(isCTFace, :);
        subs1 = currentElemLocalIdx(face2elem(isCTFace));
        subs2 = [ones(NN, 1); 2*ones(NN, 1); 3*ones(NN, 1)];
        val = B(isCTFace,:);
        for k = 1:3
            subs3 = full(currentNodelocalIdx((tFace(:, k) - 1)*CNP + subs1));
            BP = BP + accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)], val(:), [CNP, 3, nv]);
        end
    end
    
    isCSFace = isCurrentFace & ~isTriFace;
    NN = sum(isCSFace);
    if NN > 0
        sFace = face(isCSFace,:);
        subs1 = currentElemLocalIdx(face2elem(isCSFace));
        subs2 = [ones(NN, 1); 2*ones(NN, 1); 3*ones(NN, 1)];
        val = B(isCSFace,:);
        for k = 1:4
            subs3 = full(currentNodelocalIdx((sFace(:,k) - 1)*CNP + subs1));
            BP = BP + accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)], val(:), [CNP, 3, nv]);
        end
    end

    IminusP = zeros(CNP, nv, nv);
    c = 1/nv;
    ch = h(isCurrentElem);
    for i = 1:nv
        for j = 1:nv
            Xj = (node(currentElem(:, j), :) - centroid(isCurrentElem,:))./[ch, ch, ch];
            IminusP(:, i, j) = - c - dot(Xj, BP(:, :, i), 2);
            if(i == j)
                IminusP(:, i, j) = 1 + IminusP(:, i, j);
            end
        end
    end
    for i = 1:nv
        for j = 1:nv
            ii(index+1:index + CNP) = currentElem(:, i);
            jj(index+1:index + CNP) = currentElem(:, j);
            ss(index+1:index + CNP) = ch.*(dot(BP(:,:,i), BP(:,:,j),2) + dot(IminusP(:,:,i), IminusP(:,:,j),2)); 
            index = index + CNP;
        end
    end
    ft = volume(isCurrentElem).*pde.f(centroid(isCurrentElem,:))/nv;
    b = b + accumarray(currentElem(:), repmat(ft, nv, 1), [Ndof, 1]);
end

A = sparse(ii, jj, ss, Ndof, Ndof);

[fixedNode, ~, isBdNode] = findpolyboundary(face);
u = zeros(Ndof, 1);
u(fixedNode) = pde.g_D(node(fixedNode,:));
b = b - A*u;

isFreeNode = ~isBdNode;
% Set up solver
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'amg';
    end
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        tic;
        u(isFreeNode) = AD(isFreeNode,isFreeNode)\b(isFreeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'amg'
        option.solver = 'CG';
        [u(isFreeNode),info] = amg(A(isFreeNode,isFreeNode),b(isFreeNode),option);                 
end



