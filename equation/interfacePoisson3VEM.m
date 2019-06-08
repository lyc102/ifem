function [u, w, AE, AI, info] = interfacePoisson3VEM(node, elem, pde, interfaceData, option)
%% INTERFACEPOISSON3VEM Poisson equation: P1 linear element in Virtual Element Methods.
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

N = size(node, 1); Ndof = N;
NF = size(face,1); NC = size(elem, 1);

%% Deal with uncut cube elems
tic
isExteriorElem = (elemIdx == 1);
isInteriorElem = (elemIdx == -1);
exteriorElem = elem(isExteriorElem,:);
interiorElem = elem(isInteriorElem,:);

NEE = size(exteriorElem, 1);
NIE = size(interiorElem, 1);
nnzp = NEE*64;
iip = zeros(nnzp, 1);
jjp = zeros(nnzp, 1);
ssp = zeros(nnzp, 1);

nnzm = NIE*64;
iim = zeros(nnzm, 1);
jjm = zeros(nnzm, 1);
ssm = zeros(nnzm, 1);

idxp = 0;
idxm = 0;
cpp = (node(exteriorElem(:, 1),:) + node(exteriorElem(:, 7),:))/2.0;
Kp = pde.dplus(cpp);
cpm = (node(interiorElem(:, 1),:) + node(interiorElem(:, 7),:))/2.0;
Km = pde.dminus(cpm);

[AA, ~, vol] = getonecubematrix(exteriorElem(1,:));

b = zeros(Ndof, 1);
for i = 1:8
    for j = 1:8
        iip(idxp+1:idxp+NEE) = exteriorElem(:, i);
        jjp(idxp+1:idxp+NEE) = exteriorElem(:, j);
        ssp(idxp+1:idxp+NEE) = Kp*AA(i,j);
        idxp = idxp+NEE;
        
        iim(idxm+1:idxm+NIE) = interiorElem(:, i);
        jjm(idxm+1:idxm+NIE) = interiorElem(:, j);
        ssm(idxm+1:idxm+NIE) = Km*AA(i,j);
        idxm = idxm+NIE;
    end
end
ftp = vol*pde.f(cpp)/8;
ftm = vol*pde.f(cpm)/8;
for i = 1:8
    b = b + accumarray([exteriorElem(:,i);interiorElem(:,i)], [ftp;ftm], [Ndof, 1]);   
end
AE = sparse(iip, jjp, ssp, Ndof, Ndof);
AI = sparse(iim, jjm, ssm, Ndof, Ndof);
clear nnzp nnzm iim jjm ssm iip jjp ssp ftp ftm cpp Kp cpm Km AA vol



%% Deal with cutted elems

NP = max(face2elem);
isCutPoly = false(NP, 1);
isCutPoly(face2elem) = true;
cutPolyIdx = zeros(NP, 1);
isExtPoly = false(NP,1);
isExtPoly(NC+1:NP) = true;

NP = sum(isCutPoly);
cutPolyIdx(isCutPoly) = 1:NP;
face2elem = cutPolyIdx(face2elem);
isExtPoly = isExtPoly(isCutPoly);
clear isCutPoly cutPolyIdx

isTriFace = (face(:, 4) == 0);
poly2node = sparse(face2elem(isTriFace)*ones(1,3), face(isTriFace, 1:3), 1, NP,N)...
       + sparse(face2elem(~isTriFace)*ones(1,4), face(~isTriFace, 1:4), 1, NP,N);
poly2node = (poly2node > 0);
NV = poly2node*ones(N, 1);
centroid = poly2node*node./[NV, NV, NV];

coefficient = 1/3*ones(NF, 1);
coefficient(isTriFace) = 1/6;

normal = facenormal(node, face);
volume = accumarray(face2elem, coefficient.*dot(node(face(:,1),:),normal,2));
h = volume.^(1/3);

B = 1/6*normal;
B(~isTriFace, :) = 3/2*B(~isTriFace, :);
K = pde.d(centroid);

AE1 = assemblevemstiffmatrix(isExtPoly);
AI1 = assemblevemstiffmatrix(~isExtPoly);
AE = AE + AE1;
AI = AI + AI1;
A = AE + AI;
assembleTime = toc;

info.assembleTime = assembleTime;

isInterfaceNode = (vSign==0);
w = zeros(N,1);
w(isInterfaceNode) = pde.exactw(node(isInterfaceNode,:));

[fixedNode, ~, isBdNode] = findboundary3(elem);
isFreeNode = true(Ndof,1);
isFreeNode(fixedNode) = false;
u = zeros(Ndof, 1);
u(fixedNode) = pde.g_D(node(fixedNode,:));
b = b - A*u + AI*w;
b(fixedNode) = u(fixedNode);

[~, area] = facenormal(node, interface);
[lambda, weight] = quadpts(3);
nQuad = size(lambda, 1);
ge = zeros(size(interface, 1), 3);
for p = 1:nQuad
    pxy = lambda(p,1)*node(interface(:,1),:) ...
        + lambda(p,2)*node(interface(:,2),:) ...
        + lambda(p,3)*node(interface(:,3),:);
    q = pde.exactq(pxy);
    for k = 1:3
        ge(:,k) = ge(:,k) + q*lambda(p,k)*weight(p);
    end
end
ge = ge.*[area, area, area];
b = b - accumarray(interface(:), ge(:), [Ndof, 1]);

if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
switch solver
    case 'direct'
        tic;
        u(isFreeNode) = AD(isFreeNode,isFreeNode)\b(isFreeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        option.x0 = u;
        option.solver = 'CG';
        [u,info] = mg(AD,b,elem,option);
    case 'amg'
        option.solver = 'CG';
        [u(isFreeNode),info] = amg(A(isFreeNode,isFreeNode),b(isFreeNode)); 
end
info.assembleTime = assembleTime;


    function AA = assemblevemstiffmatrix(isCurrentPoly)
        
        NVN = NV(isCurrentPoly);
        if isempty(NVN)
            AA = sparse(Ndof, Ndof);
            return;
        end     
        nnz = sum(NVN.^2);
        ii = zeros(nnz, 1);
        jj = zeros(nnz, 1);
        ss = zeros(nnz, 1);
        index = 0;
        unv = unique(NVN);
        for kk = 1:length(unv)
            nv = unv(kk);
            isCurrentFace = (NV(face2elem) == nv & isCurrentPoly(face2elem));
            isCurrentElem = false(NP, 1);
            isCurrentElem(face2elem(isCurrentFace)) = true;
            % Current elem to node matrix
            currentPoly2node = poly2node(isCurrentElem, :);
            CNP = size(currentPoly2node, 1);
            currentPolyLocalIdx = zeros(NP, 1);
            currentPolyLocalIdx(isCurrentElem) = 1:CNP;
            
            
            [I, ~] = find(currentPoly2node');
            currentElem = reshape(I, nv, [])';
            localIdx = sparse(repmat((1:CNP)', 1, nv), currentElem, ones(CNP, 1)*(1:nv), CNP, N);
            clear currentPoly2node
            % Deal with current triangle face cases
            isCTFace = isCurrentFace & isTriFace;
            tFace = face(isCTFace, :);
            NN = sum(isCTFace);
            subs1 = currentPolyLocalIdx(face2elem(isCTFace));
            subs2 = [ones(NN, 1); 2*ones(NN, 1); 3*ones(NN, 1)];
            val = B(isCTFace,:);
            BP = zeros(CNP,3,nv);
            for m = 1:3
                subs3 = full(localIdx((tFace(:, m) - 1)*CNP + subs1));
                BP = BP + accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)], val(:), [CNP, 3, nv]);
            end
 
            isCSFace = isCurrentFace & ~isTriFace;
            sFace = face(isCSFace,:);
            NN = size(sFace,1);
            if NN > 0
                subs1 = currentPolyLocalIdx(face2elem(isCSFace));
                subs2 = [ones(NN, 1); 2*ones(NN, 1); 3*ones(NN, 1)];
                val = B(isCSFace,:);
                for m = 1:4
                    subs3 = full(localIdx((sFace(:, m) - 1)*CNP + subs1));
                    BP = BP + accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)], val(:), [CNP, 3, nv]);
                end
            end
            clear isCurrentFace currentPolyLocalIdx localIdx
            clear subs1 subs2 subs3
            
            IminusP = zeros(CNP, nv, nv);
            for n = 1:nv
                for m = 1:nv
                    Xn = node(currentElem(:, n), :) - centroid(isCurrentElem,:);
                    IminusP(:, n, m) = - 1/nv - dot(Xn, BP(:, :, m), 2)./volume(isCurrentElem);
                    if(n == m)
                        IminusP(:, n, m) = 1 + IminusP(:, n, m);
                    end
                end
            end
            for n = 1:nv
                for m = 1:nv
                    ii(index+1:index + CNP) = currentElem(:, n);
                    jj(index+1:index + CNP) = currentElem(:, m);
                    ss(index+1:index + CNP) = K(isCurrentElem).*(dot(BP(:,:,n), BP(:,:,m),2)./volume(isCurrentElem) ...
                        + h(isCurrentElem).*dot(IminusP(:,:,n), IminusP(:,:,m),2));
                    index = index + CNP;
                end
            end
            clear IminusP BP
            ft = volume(isCurrentElem).*pde.f(centroid(isCurrentElem,:))/nv;
            b = b + accumarray(currentElem(:), repmat(ft, nv, 1), [Ndof, 1]);
            clear ft isCurrentElem currentElem 
        end  
        AA = sparse(ii, jj, ss, Ndof, Ndof);
    end
  

    function [AA, hh, vol] = getonecubematrix(cubeElem)
        localIdxFace = [1, 4, 3, 2
                        6, 7, 8, 5
                        2, 3, 7, 6
                        1, 5, 8, 4
                        1, 2, 6, 5
                        4, 8, 7, 3];
        cubeFace = cubeElem(localIdxFace);
        cp = (node(cubeElem(1), :) + node(cubeElem(7),:))/2.0;
        fn = facenormal(node,cubeFace);
        vol = sum(dot(node(cubeFace(:,1),:),fn,2)/3);
        hh = vol^(1/3);
        BB = 0.25*fn;       
        BN = zeros(8, 3);
        subs2 = [ones(6,1);2*ones(6,1);3*ones(6,1)];
        for pp = 1:4
           subs1 = localIdxFace(:,pp);
           BN = BN + accumarray([repmat(subs1,3,1),subs2], BB(:),[8,3]);
        end
        D = node(cubeElem,:) - repmat(cp, 8, 1);
        IP = eye(8)- 1/8 - D*BN'/vol;
        AA = BN*BN'/vol + hh*(IP'*IP);
    end
end
