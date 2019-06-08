function [u,p,info,Ai,Bi,invMt,Pro_u,Pro_p] = mgstokesRT0(A,B,f,g,u,p,node,elem,ufreeDof,option,varargin) 
%% MGSTOKES 
%
% Created by Ming Wang and Jie Zhou based on discussion with Long Chen.
% Improvement by Long Chen: improve efficiency, include smoother and add
% bisection case and BDM1B.

% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

tic;

% Assign default values to unspecified parameters
if ~exist('option','var'), 
    option = []; 
end
option = mgoptions(option,length(f)+length(g));    % parameters
tol = option.tol; maxIt = option.solvermaxit; 
smoothingstep = option.smoothingstep; solver = option.solver;
printlevel = option.printlevel; setupflag = option.setupflag;

if setupflag == true
%% Reconstruced hierarchical meshes and transfer operators
N0 = 8;
N = size(node,1); NT = size(elem,1);
level = 20;
NL(level+1) = N; % now NL(1:level) = 0;
nodei = cell(level,1);
elemi = cell(level,1);
Pro_u = cell(level,1);
Pro_p = cell(level,1);
nodei{level} = node;
elemi{level} = elem;
for j = level: -1 : 2
    switch option.refType 
        case 'red'
            [elemi{j-1},newHB] = uniformcoarsenred(elemi{j}); % coasen red refinement
            if ~isempty(newHB)
                Pro_u{j-1} = transferedgered(elemi{j-1},elemi{j}); % transfer operator of u
                Pro_p{j-1} = repmat(speye(size(elemi{j-1},1)),4,1);   % transfer operator of p
            end
        case 'bisect'
            % merge two coarsen of bisection grids s.t. the ratio is 1/4.
            [tempelem,newHB,tree] = uniformcoarsen(elemi{j}); % coarse bisection
            if ~isempty(newHB)            
                % first coarsen
                tempPro_u = transferedgecoarsen(tempelem,elemi{j},tree);
                tempPro_p = transferelem(tempelem,elemi{j},tree);
                % second coarsen
                [elemi{j-1},newHB,tree] = uniformcoarsen(tempelem); % coarse bisection
                Pro_u{j-1} = transferedgecoarsen(elemi{j-1},tempelem,tree);
                Pro_u{j-1} = tempPro_u*Pro_u{j-1};
                Pro_p{j-1} = transferelem(elemi{j-1},tempelem,tree);
                Pro_p{j-1} = tempPro_p*Pro_p{j-1};
            end
    end    
    if (isempty(newHB)) || (size(elemi{j-1},1)< 2*N0) 
    % no nodes are removed or it reaches the coarsest level
        NL = NL(j:end);       
        break; 
    end
    NL(j) = NL(j+1) - size(newHB,1); % update NL(k)
    nodei{j-1} = nodei{j};           % node in i-th level
    nodei{j-1}(newHB(:,1),:) = [];   % remove fine nodes
end
level = length(NL)-1;    % actual level
nodei = nodei(end-level+1: end);
elemi = elemi(end-level+1: end);
Pro_u = Pro_u(end-level+1: end);
Pro_p = Pro_p(end-level+1: end);

%% Matrices in each level
Ai = cell(level,1);
Ai{level} = A;
Bi = cell(level,1);
Bi{level} = B; 
invMt = cell(level,1);
area = simplexvolume(node,elem);
NTf = size(elem,1); 
invMt{level} = spdiags(1./area,0,NTf,NTf);
NE = size(Pro_u{level-1},1);
ufreeDofj = ufreeDof(ufreeDof <= NE);
for j = level:-1:2   
    % get mesh and data in the j-1 level
    elem = sortelem(elemi{j-1});
    node = nodei{j-1};
    [elem2edge,edge] = dofedge(elem);
    [Dlambda,area,elemSign] = gradbasis(node,elem);
    NTj = size(elem,1); 
    
    % assemble matrix using the mesh information.
    % Mv: Lumped mass matrix for vertex: P1 element
    Mv = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
    Mv((Mv == 0)) = [];
    Nj = length(Mv); 
    invMv = spdiags(1./Mv,0,Nj,Nj);
    Mv = spdiags(Mv,0,Nj,Nj);
    % Me: Mass matrix for RT0 element
    Me = getmassmatvec(elem2edge,area,Dlambda,'RT0');
    %invMt: the inverse of Mass matrix for P0 element
    invMt{j-1} = spdiags(1./area,0,NTj,NTj);
    % B: - divergence operator
    Bj = -icdmat(double(elem2edge),elemSign*[1 -1 1]);
    % C: curl operator
    C = icdmat(double(edge),[-1 1]);
    % R: weak rot operator
    R = invMv*C'*Me;
    % Vector Laplacian
    Ai{j-1} =  Bj'*invMt{j-1}*Bj + R'*Mv*R;
    
    % find free dof in each level
    isFixedDof = false(size(Ai{j-1},1),1);
    isDirichlet = false(size(Ai{j-1},1),1);
    % add coarsening into bdFlag
    bdFlag = setboundary(node,elem,'Dirichlet');
    if ~isempty(bdFlag)       
        isDirichlet(elem2edge(bdFlag(:)==1)) = true;
        isFixedDof(isDirichlet) = true;
        ufreeDofc = find(~isFixedDof);
    end
    
    % modify transfer operators
    Pro_u{j-1} = Pro_u{j-1}(ufreeDofj,ufreeDofc);
    Ai{j-1} = Ai{j-1}(ufreeDofc,ufreeDofc);
    Bi{j-1} = Bj(:,ufreeDofc);
    ufreeDofj = ufreeDofc;
end
clear nodei elemi

%% Modify the prolongation operator for BDM1B
if length(u) == 2*NE+NT
    % add one more level from RT0 to BDM1
    NRT0 = size(Pro_u{level-1},1);
    Pro_u{level-1} = [Pro_u{level-1}; sparse(size(A,1)-NRT0,size(Ai{level-1},1))];
    Res_u{level} = Pro_u{level-1}';
end

end % end for setup

%% Only set up the transfer operators
if strcmp(solver,'NO') 
    flag = 0; itStep = 0; err = 0; time = toc;
    info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));
    return
end

%% No need of set up
if setupflag == false
    Ai = varargin{1};
    Bi = varargin{2};
   invMt = varargin{3};
   Pro_u = varargin{4};
   Pro_p = varargin{5};
   level = size(Ai,1);
end

%% Construct other necessary operators
Nu = zeros(level,1);
Nu(level) = size(A,1);                
Np = zeros(level,1);
Np(level) = size(B,1);
bigAi = cell(level,1);
bigAi{level} = [Ai{level} Bi{level}'; Bi{level} sparse(Np(level),Np(level))];
Res_u = cell(level,1);
Res_p = cell(level,1);
for j = level:-1:2   
    Res_u{j} = Pro_u{j-1}';
    Res_p{j} = Pro_p{j-1}';    
    Nu(j-1) = size(Ai{j-1},1); 
    Np(j-1) = size(Bi{j-1},1);    
    bigAi{j-1} = [Ai{j-1} Bi{j-1}'; Bi{j-1} sparse(Np(j-1),Np(j-1))];
end
Ndof = Nu + Np;        

%% Smoothers
Su = cell(level,1);
Sp = cell(level,1);
for k = 1:level
    Sp{k} = tril(Bi{k}*(Bi{k})');
    Su{k} = tril(Ai{k});
end

%% Solver
bigF = [f; g-mean(g)];
bigu = [u(ufreeDof); p];
bigr = bigF - bigAi{level}*bigu;
nb = norm(bigF);
err = zeros(maxIt,1);
err(1) = norm(bigr)/nb;
k = 1;
% condest(bigAi{J}(1:end-1,1:end-1))
if strcmp(solver,'GMRES') % GMRES solver
    restart = 50; prefunc = @(r)kcycle(r,level);
    tic;
    [bige, itStep, relRes_u] = PFGMRES(bigAi{level}, bigr, bigu, maxIt, restart, tol, prefunc,0);
    bigu = bigu + bige;
    time = toc;
    fprintf('itStep=%d, relRes_u = %10.6g, time = %4.3g s,\n', ...
             itStep, relRes_u(itStep), time);
else  % MG
    while (max(err(k)) > tol) && (k <= maxIt)
        k = k + 1;
        % one multigrid cycle
        switch (solver)
            case 'VCYCLE'
                biguerr = vcycle(bigr);
            case 'WCYCLE'
                biguerr = wcycle(bigr);
        end
        % update solution
        bigu = bigu + biguerr;
        % compute residual
        bigr = bigr - bigAi{level}*biguerr;
        % compute the relative errror
        err(k) = norm(bigr)/nb;
        if printlevel >= 2
            fprintf('#dof: %8.0u, MG %8s iter: %2.0u, err = %8.4e\n',...
                 Ndof(level), solver, k-1, err(k));
        end
    end
end
err = err(1:k);
itStep = k-1;
u = bigu(1:Nu(level)); 
p = bigu(Nu(level)+1:end);

%% Output
if k > maxIt
    flag = 1;
else
    flag = 0;
end
time = toc;
if printlevel >= 1
    fprintf('#dof: %6.0u,  #nnz: %6.0u, level: %2.0u  MG %6s iter: %2.0u,  err = %8.4e,  time = %4.2g s\n',...
             Ndof(level), nnz(bigAi{level}), level, solver, itStep, err(end), time)
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions vcycle, wcycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Vcycle MG
    function e = vcycle(r,J)
        if nargin < 2, J = level; end
        if J == 1 % exact solver in the coaresest grid
           e = zeros(size(r));
           e(1:end-1) = bigAi{J}(1:end-1,1:end-1)\r(1:end-1);
           e(Nu(J)+1:end) = e(Nu(J)+1:end)-mean(e(Nu(J)+1:end));                    
           return
        end
        ru = r(1:Nu(J)); 
        rp = r(Nu(J)+1:end);
        
        % pre-smoothing in the fine grid 
        [eu,ep] = DGSRT0(Ai{J},Bi{J},Bi{J}',zeros(Nu(J),1),zeros(Np(J),1),...
                         ru,rp,smoothingstep,Su{J},Sp{J},invMt{J});
        
        % form residual and restrict onto coarse grid
        ruc = Res_u{J}*(ru - Ai{J}*eu-(Bi{J})'*ep);
        rpc = Res_p{J}*(rp - Bi{J}*eu);          
        
        % coarse grid correction twice
        rc = [ruc;rpc];
        ec = vcycle(rc,J-1);

        % prolongate coarse grid correction to the fine grid
        tempeu = Pro_u{J-1}*ec(1:Nu(J-1));
        tempep = Pro_p{J-1}*ec(Nu(J-1)+1:end);    
        eu = eu + tempeu; 
        ep = ep + tempep;

        % post-smoothing in the fine grid
        [eu,ep] = DGSRT0(Ai{J},Bi{J},Bi{J}',eu,ep,ru,rp,smoothingstep,...
                         Su{J},Sp{J},invMt{J});
        e = [eu;ep];
    end

%% Wcycle MG
    function e = wcycle(r,J)
        if nargin < 2, J = level; end
        if J == 1 % exact solver in the coaresest grid
           e = zeros(size(r));
           e(1:end-1) = bigAi{J}(1:end-1,1:end-1)\r(1:end-1);
           e(Nu(J)+1:end) = e(Nu(J)+1:end)-mean(e(Nu(J)+1:end));                    
           return
        end
        ru = r(1:Nu(J)); 
        rp = r(Nu(J)+1:end);
        
        % pre-smoothing in the fine grid 
        [eu,ep] = DGSRT0(Ai{J},Bi{J},Bi{J}',zeros(Nu(J),1),zeros(Np(J),1),...
                         ru,rp,smoothingstep,Su{J},Sp{J},invMt{J});
        
        % form residual and restrict onto coarse grid
        ruc = Res_u{J}*(ru - Ai{J}*eu-(Bi{J})'*ep);
        rpc = Res_p{J}*(rp - Bi{J}*eu);          
        
        % coarse grid correction twice
        rc = [ruc;rpc];
        ec = wcycle(rc,J-1);
        ec = ec + wcycle(rc-bigAi{J-1}*ec,J-1);

        % prolongate coarse grid correction to the fine grid
        tempeu = Pro_u{J-1}*ec(1:Nu(J-1));
        tempep = Pro_p{J-1}*ec(Nu(J-1)+1:end);    
        eu = eu + tempeu; 
        ep = ep + tempep;

        % post-smoothing in the fine grid
        [eu,ep] = DGSRT0(Ai{J},Bi{J},Bi{J}',eu,ep,ru,rp,smoothingstep,...
                         Su{J},Sp{J},invMt{J});
        e = [eu;ep];
    end
 
%% DGS for RT0
    function [u,p] = DGSRT0(A,B,Bt,u,p,f,g,itStep,Su,Sp,invMp)
        for s = 1: itStep
            % Step 1: relax Momentum eqns by G-S
            for i = 1:2
                u = u + Su\(f-Bt*p-A*u);
            end
            % Step 2: relax Continuity eqns by G-S
            rp = g - B*u;
%             dq = Sp\rp;
            % Step 3: update u and p
%             u = u + Bt*dq;
            p = p - 1.2*invMp*rp;
%             p = p - Sp\(B*(A*(Bt*dq)));
        end
    end
end