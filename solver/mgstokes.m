function [u,p,info] = mgstokes(A,B,f,g,u,p,elem,freeDof,option) 
%% MGSTOKES: multigrid methods for Stokes equations
%
% Created by Ming Wang and Jie Zhou on Mar,3,2013 for P2 elements. Add
% other elements by Jie and Long on Mar 24,2013.
%
% Reference: 
%
% M. Wang and L. Chen. Multigrid Methods for the Stokes equations using
% Distributive Gauss-Seidel Relaxations based on the Least Squares
% Commutator. Journal of Scientific Computing. 56(2): 409-431, 2013.
%
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

t = cputime;

%% Size of systems
N = max(elem(:));
Ndof = length(f) + length(g);

%% Options
% Assign default values to unspecified parameters
if ~exist('option','var')
    option = []; 
end
option = mgoptions(option,Ndof);    % parameters
tol = option.tol; maxIt = option.solvermaxit; 
printlevel = option.printlevel; 
solver = option.solver;

%% Mesh on each level
N0 = 10;
HB = zeros(N,3);
level = 20;
NL(level+1) = N; % now NL(1:level) = 0;
elemi = cell(level,1);
elemi{level} = elem;
for j = level: -1 : 2
    [elemi{j-1},newHB] = uniformcoarsenred(elemi{j});  % try coasen red refinement
    if (isempty(newHB)) || (size(elemi{j-1},1)< 2*N0) 
    % no nodes are removed or it reaches the coarsest level
        NL = NL(j:end);       
        break; 
    end
    NL(j) = NL(j+1) - size(newHB,1); % update NL(k)
    HB(NL(j)+1:NL(j+1),1:3) = newHB(:,1:3);
end
level = length(NL)-1;    % actual level
elemi = elemi(end-level+1:end);

%% Transfer operators between levels
Pro_u = cell(level,1);
Res_u = cell(level,1);
Pro_p = cell(level,1);
Res_p = cell(level,1);

switch(upper(option.elemType(1:end-2))) % transfer operator of velocity 
    case 'P2'
        P2FreeDof = freeDof(1:end/2);   % again only one component
        for i = level:-1:2 
            [ProU0,P2FreeDof] = transferP2red(elemi{i-1},elemi{i},P2FreeDof);
            Pro_u{i-1} = ProU0;
            Res_u{i}   = ProU0';
        end         
    case 'CR'
        CRfreeDof = freeDof(1:end/2);
        for i=level:-1:2 
            [ProU0,CRfreeDof] = transferCRred(elemi{i-1},elemi{i},CRfreeDof);
            Pro_u{i-1} = ProU0;
            Res_u{i}   = ProU0';
        end     
    case 'P1B'    
        for i=level:-1:2 
            Pro_u{i-1} = speye(length(freeDof),length(freeEdge));
            Res_u{i} = Pro_u{i-1}';
        end
    case 'ISOP2'
        P2FreeDof = freeDof(1:end/2);  % again only one component
        for i = level:-1:2 
            [ProU0,P2FreeDof] = nodeisoP2TansferUniform(elemi{i-1},elemi{i},P2FreeDof);
            Pro_u{i-1} = ProU0;
            Res_u{i}   = ProU0';
        end                 
end
switch upper(option.elemType(end-1:end))    
    case  'P0'
        for k = level:-1:2 
            Np(k-1)    = size(elemi{k-1},1);
            Pro_p{k-1} = repmat(speye(Np(k-1)),4,1);
            Res_p{k} = Pro_p{k-1}';
        end   
    case  'P1'
       [Pro_p,Res_p] = transferoperator(HB,NL); 
end

%% Matrices in each level
% size of u and p
Nu = zeros(level,1); 
Np = zeros(level,1);
Nu(level) = size(A,1)/2; % one component of velocity
Np(level) = size(B,1);

Ai = cell(level,1);
Bi = cell(level,1);
bigAi = cell(level,1);
Ai{level} = A; 
Bi{level} = B;
bigAi{level} = [Ai{level} Bi{level}'; Bi{level} sparse(Np(level),Np(level))];

for j = level:-1:2
    Ai{j-1} = blkdiag(Res_u{j},Res_u{j})*Ai{j}*blkdiag(Pro_u{j-1},Pro_u{j-1});  % Ac = Res*Af*Pro
    Bi{j-1} = Res_p{j}*Bi{j}*blkdiag(Pro_u{j-1},Pro_u{j-1});
    Nu(j-1) = size(Ai{j-1},1)/2; 
    Np(j-1) = size(Bi{j-1},1);
    bigAi{j-1} = [Ai{j-1} Bi{j-1}'; Bi{j-1} sparse(Np(j-1),Np(j-1))];
end
Ndof = 2*Nu+Np;

%% matrix used for smoothing
if isfield(option,'smootherType')
    smoother = option.smootherType;
else 
    smoother ='LSCDGS';
end
auxMat = cell(level,1);
switch(smoother)
    case 'LSCDGS'
        for k = 2:level
            Bt = (Bi{k})';
            BBt = Bi{k}*Bt;
            BABt = Bi{k}*Ai{k}*Bt;
            Su = tril(Ai{k});
            Sp = tril(BBt);
            Spt = triu(BBt);
            DSp = diag(BBt);
            auxMat{k} = struct('Bt',Bt,'BBt',BBt,'BABt',BABt,'Su',Su,...
                               'Spt',Spt,'Sp',Sp,'DSp',DSp);
        end
    case 'IUzawa'
        for k = 2:level
            Bt = (Bi{k})';
            DA = 2*diag(Ai{k});
            invDA = spdiags(1./DA,0,2*Nu(k),2*Nu(k));
            BinvDABt = Bi{k}*invDA*Bt;
            auxMat{k} = struct('Bt',Bt,'DA',DA,'BinvDABt',BinvDABt);
        end
end
Ai_IU = cell(level,1);
Si_IU = cell(level,1);
SSi_IU = cell(level,1);
Res_IU = cell(level,1);
Pro_IU = cell(level,1);
if ~isfield(option,'smootherbarSp')
    option.smootherbarSp ='SGS';
end
if strcmp(option.smootherbarSp,'VCYCLE')
%     option.smootherbarSp = 'VCYCLE';
    innerMGoption.solver = 'NO';
    innerMGoption.N0 = 10;
    for j = 2:level
        [tempvar,tempvar,Ai_IU{j},Si_IU{j},SSi_IU{j},Res_IU{j},Pro_IU{j}] = ...
                mg(auxMat{j}.BBt,ones(Np(j),1),elemi{j},innerMGoption); 
    end
end

%% Solver
%  initial set up
bigF = [f; g-mean(g)];
bigu = [u; p];
bigr = bigF - bigAi{level}*bigu;
nb = norm(bigF);
err = zeros(maxIt,1);
err(1) = norm(bigr)/nb;
k = 1;

while (max(err(k)) > tol) && (k <= maxIt)
    k = k + 1;
    switch (solver)
        case 'VCYCLE'
            bigerru = vcycle(bigr);
        case 'WCYCLE'
            bigerru = wcycle(bigr);
    end
    bigu = bigu+bigerru;
    bigr = bigr - bigAi{level}*bigerru;
    % compute the relative error
    err(k) = norm(bigr)/nb;
    if printlevel >= 2    
        fprintf('#dof: %8.0u, MG %8s iter: %2.0u, err = %8.4e\n',...
            Ndof(level), solver,k-1, err(k));
    end
end
err = err(1:k);
itStep = k-1;
u = bigu(1:2*Nu(level));  
p = bigu(2*Nu(level)+1:end);

%% Output
if k > maxIt
    flag = 1;
else
    flag = 0;
end
time = cputime-t;
if printlevel >= 1
    fprintf('#dof: %6.0u,  #nnz: %6.0u, level: %2.0u  MG %6s iter: %2.0u,  err = %8.4e,  time = %4.2g s\n',...
             Ndof(level), nnz(bigAi{level}), level, solver, itStep, err(end), time);    
end
if (flag == 1) && (printlevel>0)
   fprintf('NOTE: the iterative method does not converge! \n');    
end
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions vcycle, wcycle, smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% subfunctions v-cycle
     function e = vcycle(r,J)
        if nargin < 2, J = level; end
        if J == 1 % solver in the coaresest grid
            e = zeros(size(r));
            e(1:end-1) = bigAi{J}(1:end-1,1:end-1)\r(1:end-1);
            e(2*Nu(J)+1:end) = e(2*Nu(J)+1:end)-mean(e(2*Nu(J)+1:end));
            return
        end
        ru = r(1:2*Nu(J)); rp = r(2*Nu(J)+1:end);

        % pre-smoothing
        [eu,ep] = smoothing(zeros(2*Nu(J),1),zeros(Np(J),1),ru,rp,J);
        
        % form residual and restrict onto coarse grid
        rru = ru-Ai{J}*eu-(Bi{J})'*ep;
        rrp = rp-Bi{J}*eu;
        ruc = reshape(Res_u{J}*reshape(rru,Nu(J),2),2*Nu(J-1),1);
        rpc = Res_p{J}*rrp;
        
        % coarse grid correction
        rc = [ruc;rpc];        
        ec = vcycle(rc,J-1);         
        
        % correction on the fine grid
        tempeu = reshape(Pro_u{J-1}*reshape(ec(1:2*Nu(J-1)),Nu(J-1),2),2*Nu(J),1);
        tempep = Pro_p{J-1}*ec(2*Nu(J-1)+1:end);
        eu = tempeu + eu; 
        ep = tempep + ep;
        
        % post-smoothing
        [eu,ep] = smoothing(eu,ep,ru,rp,J);
        e = [eu;ep];

     end
 
 %% Wcycle MG
     function e = wcycle(r,J)
        if nargin < 2, J = level; end
        if J == 1 % solver in the coaresest grid
            e = zeros(size(r));
            e(1:end-1) = bigAi{J}(1:end-1,1:end-1)\r(1:end-1);
            e(2*Nu(J)+1:end) = e(2*Nu(J)+1:end)-mean(e(2*Nu(J)+1:end));
            return
        end
        ru = r(1:2*Nu(J)); 
        rp = r(2*Nu(J)+1:end);
%       
        % pre-smoothing in the fine grid 
        [eu,ep] = smoothing(zeros(2*Nu(J),1),zeros(Np(J),1),ru,rp,J);
        
        % form residual and restrict onto coarse grid
        rru = ru - Ai{J}*eu-(Bi{J})'*ep;
        rrp = rp - Bi{J}*eu;
        ruc = reshape(Res_u{J}*reshape(rru,Nu(J),2),2*Nu(J-1),1);
        rpc = Res_p{J}*rrp;

        % coarse grid correction
        rc = [ruc;rpc];
        ec = wcycle(rc,J-1);
        % once more for w-cycle
        ec = ec+ wcycle(rc - bigAi{J-1}*ec,J-1);
        
        % correction on the fine grid
        tempeu = reshape(Pro_u{J-1}*reshape(ec(1:2*Nu(J-1)),Nu(J-1),2),2*Nu(J),1);
        tempep = Pro_p{J-1}*ec(2*Nu(J-1)+1:end);
        eu = eu + tempeu; 
        ep = ep + tempep;
        
        % post-smoothing in the fine grid
        [eu,ep] = smoothing(eu,ep,ru,rp,J);
        e = [eu;ep];
     end
 
 %% Smoothing
    function [u,p] = smoothing(u,p,f,g,J)
        if strcmp(smoother,'LSCDGS')
            [u,p] = StokesLSCDGS(u,p,f,g,Ai{J},Bi{J},auxMat{J},elemi{J},option,...
                                 Ai_IU{J},Si_IU{J},SSi_IU{J},Res_IU{J},Pro_IU{J});
        elseif strcmp(smoother,'BDDGS')
            % e = StokesBiDGStri(Ai{J}, Bi{J}, e, r, smoothStep, invDM{J},area{J});
        elseif strcmp(smoother,'iUzawa')
            [u,p] = StokesIUzawa(u,p,f,g,Ai{J},Bi{J},auxMat{J},elemi{J},option,...
                                 Ai_IU{J},Si_IU{J},SSi_IU{J},Res_IU{J},Pro_IU{J});
        end
    end
end