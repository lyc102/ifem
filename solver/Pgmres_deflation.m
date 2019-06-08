function [x,flag,relres,iter,resvec] = Pgmres_deflation(A,b,restart,tol,maxit,M1,M2,x,varargin)
%GMRES   Generalized Minimum Residual Method.
%   X = GMRES(A,B) attempts to solve the system of linear equations A*X = B
%   for X.  The N-by-N coefficient matrix A must be square and the right
%   hand side column vector B must have length N. This uses the unrestarted
%   method with MIN(N,10) total iterations.
%
%   X = GMRES(AFUN,B) accepts a function handle AFUN instead of the matrix
%   A. AFUN(X) accepts a vector input X and returns the matrix-vector
%   product A*X. In all of the following syntaxes, you can replace A by
%   AFUN.
%
%   X = GMRES(A,B,RESTART) restarts the method every RESTART iterations.
%   If RESTART is N or [] then GMRES uses the unrestarted method as above.
%
%   X = GMRES(A,B,RESTART,TOL) specifies the tolerance of the method.  If
%   TOL is [] then GMRES uses the default, 1e-6.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT) specifies the maximum number of outer
%   iterations. Note: the total number of iterations is RESTART*MAXIT. If
%   MAXIT is [] then GMRES uses the default, MIN(N/RESTART,10). If RESTART
%   is N or [] then the total number of iterations is MAXIT.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M) and
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M1,M2) use preconditioner M or M=M1*M2
%   and effectively solve the system inv(M)*A*X = inv(M)*B for X. If M is
%   [] then a preconditioner is not applied.  M may be a function handle
%   returning M\X.
%
%   X = GMRES(A,B,RESTART,TOL,MAXIT,M1,M2,X0) specifies the first initial
%   guess. If X0 is [] then GMRES uses the default, an all zero vector.
%
%   [X,FLAG] = GMRES(A,B,...) also returns a convergence FLAG:
%    0 GMRES converged to the desired tolerance TOL within MAXIT iterations.
%    1 GMRES iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 GMRES stagnated (two consecutive iterates were the same).
%
%   [X,FLAG,RELRES] = GMRES(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL. Note with
%   preconditioners M1,M2, the residual is NORM(M2\(M1\(B-A*X))).
%
%   [X,FLAG,RELRES,ITER] = GMRES(A,B,...) also returns both the outer and
%   inner iteration numbers at which X was computed: 0 <= ITER(1) <= MAXIT
%   and 0 <= ITER(2) <= RESTART.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = GMRES(A,B,...) also returns a vector of
%   the residual norms at each inner iteration, including NORM(B-A*X0).
%   Note with preconditioners M1,M2, the residual is NORM(M2\(M1\(B-A*X))).
%
%   Example:
%      n = 21; A = gallery('wilk',n);  b = sum(A,2);
%      tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%      x = gmres(A,b,10,tol,maxit,M);
%   Or, use this matrix-vector product function
%      %-----------------------------------------------------------------%
%      function y = afun(x,n)
%      y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%      %-----------------------------------------------------------------%
%   and this preconditioner backsolve function
%      %------------------------------------------%
%      function y = mfun(r,n)
%      y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%      %------------------------------------------%
%   as inputs to GMRES:
%      x1 = gmres(@(x)afun(x,n),b,10,tol,maxit,@(x)mfun(x,n));
%
%   Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%      float: double
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, LSQR, MINRES, PCG, QMR, SYMMLQ,
%   TFQMR, ILU, FUNCTION_HANDLE.

%   References
%   H.F. Walker, "Implementation of the GMRES Method Using Householder
%   Transformations", SIAM J. Sci. Comp. Vol 9. No 1. January 1988.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.21.4.15 $ $Date: 2011/05/17 02:33:07 $
global ZH AH
if (nargin < 2)
    error(message('MATLAB:gmres:NumInputs'));
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);
if strcmp(atype,'matrix')
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        error(message('MATLAB:gmres:SquareMatrix'));
    end
    if ~isequal(size(b),[m,1])
        error(message('MATLAB:gmres:VectorSize', m));
    end
else
    m = size(b,1);
    n = m;
    if ~iscolumn(b)
        error(message('MATLAB:gmres:Vector'));
    end
end

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(restart) || (restart == n)
    restarted = false;
else
    restarted = true;
end
if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:gmres:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:gmres:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 5) || isempty(maxit)
    if restarted
        maxit = min(ceil(n/restart),10);
    else
        maxit = min(n,10);
    end
end

if restarted
    outer = maxit;
    if restart > n
        warning(message('MATLAB:gmres:tooManyInnerItsRestart',restart, n));
        restart = n;
    end
    inner = restart;
else
    outer = 1;
    if maxit > n
        warning(message('MATLAB:gmres:tooManyInnerItsMaxit',maxit, n));
        maxit = n;
    end
    inner = maxit;
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                   % Norm of rhs vector, b
if (n2b == 0)                    % if    rhs vector is all zeros
    x = zeros(n,1);              % then  solution is all zeros
    flag = 0;                    % a valid solution has been obtained
    relres = 0;                  % the relative residual is actually 0/0
    iter = [0 0];                % no iterations need be performed
    resvec = 0;                  % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('gmres',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if ((nargin >= 6) && ~isempty(M1))
    existM1 = 1;
    [m1type,m1fun,m1fcnstr] = iterchk(M1);  %Check the size of Precondioner.
    if strcmp(m1type,'matrix')
        if ~isequal(size(M1),[m,m])
            error(message('MATLAB:gmres:PreConditioner1Size', m));
        end
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 7) && ~isempty(M2))
    existM2 = 1;
    [m2type,m2fun,m2fcnstr] = iterchk(M2);
    if strcmp(m2type,'matrix')
        if ~isequal(size(M2),[m,m])
            error(message('MATLAB:gmres:PreConditioner2Size', m));
        end
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

if ((nargin >= 8) && ~isempty(x))
    if ~isequal(size(x),[n,1])
        error(message('MATLAB:gmres:XoSize', n));
    end
else
    x = zeros(n,1);
end

if ((nargin > 8) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:gmres:TooManyInputs'));
end

% Set up for the method
flag = 1;
xmin = x;                        % Iterate which has minimal residual so far
imin = 0;                        % "Outer" iteration at which xmin was computed
jmin = 0;                        % "Inner" iteration at which xmin was computed
tolb = tol * n2b;                % Relative tolerance
evalxm = 0;
stag = 0;
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;
minupdated = 0;

x0iszero = (norm(x) == 0);
r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
normr = norm(r);                 % Norm of initial residual
if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = [0 0];
    resvec = normr;
    if (nargout < 2)
        itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
    end
    return
end
minv_b = b;

if existM1
    
%    r = P*r;
    r = r-A*(ZH*(AH\(ZH'*r)));
    r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
    if ~x0iszero
        
        b0 = b-A*(ZH*(AH\(ZH'*b)));
        minv_b = iterapp('mldivide',m1fun,m1type,m1fcnstr,b0,varargin{:});
    else
        minv_b = r;
    end
    if ~all(isfinite(r)) || ~all(isfinite(minv_b))
        flag = 2;
        x = xmin;
        relres = normr / n2b;
        iter = [0 0];
        resvec = normr;
        return
    end
end

if existM2
    r = iterapp('mldivide',m2fun,m2type,m2fcnstr,r,varargin{:});
    if ~x0iszero
        minv_b = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_b,varargin{:});
    else
        minv_b = r;
    end
    if ~all(isfinite(r)) || ~all(isfinite(minv_b))
        flag = 2;
        x = xmin;
        relres = normr / n2b;
        iter = [0 0];
        resvec = normr;
        return
    end
end

normr = norm(r);                 % norm of the preconditioned residual
n2minv_b = norm(minv_b);         % norm of the preconditioned rhs
clear minv_b;
tolb = tol * n2minv_b;
if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2minv_b;
    iter = [0 0];
    resvec = n2minv_b;
    if (nargout < 2)
        itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
    end
    return
end

resvec = zeros(inner*outer+1,1);  % Preallocate vector for norm of residuals
resvec(1) = normr;                % resvec(1) = norm(b-A*x0)
normrmin = normr;                 % Norm of residual from xmin

%  Preallocate J to hold the Given's rotation constants.
J = sparse(2,inner);

U = sparse(n,inner);
R = sparse(inner,inner);
w = sparse(inner+1,1);

for outiter = 1 : outer
    %  Construct u for Householder reflector.
    %  u = r + sign(r(1))*||r||*e1
    u = r;
    normr = norm(r);
    beta = scalarsign(r(1))*normr;
    u(1) = u(1) + beta;
    u = u / norm(u);
    
    U(:,1) = u;
    
    %  Apply Householder projection to r.
    %  w = r - 2*u*u'*r;
    w(1) = -beta;
    
    for initer = 1 : inner
        %  Form P1*P2*P3...Pj*ej.
        %  v = Pj*ej = ej - 2*u*u'*ej
        v = -2*(u(initer)')*u;
        v(initer) = v(initer) + 1;
        %  v = P1*P2*...Pjm1*(Pj*ej)
        for k = (initer-1):-1:1
            v = v - U(:,k)*(2*(U(:,k)'*v));
        end
        %  Explicitly normalize v to reduce the effects of round-off.
        v = v/norm(v);
        
        %  Apply A to v.
        v = iterapp('mtimes',afun,atype,afcnstr,v,varargin{:});
        %  Apply Preconditioner.
        if existM1
            
            v = v-A*(ZH*(AH\(ZH'*v)));
            v = iterapp('mldivide',m1fun,m1type,m1fcnstr,v,varargin{:});
            %v = P*v;
            if ~all(isfinite(v))
                flag = 2;
                break
            end
        end
        
        if existM2
            v = iterapp('mldivide',m2fun,m2type,m2fcnstr,v,varargin{:});
            if ~all(isfinite(v))
                flag = 2;
                break
            end
        end
        %  Form Pj*Pj-1*...P1*Av.
        for k = 1:initer
            v = v - U(:,k)*(2*(U(:,k)'*v));
        end
        
        %  Determine Pj+1.
        if (initer ~= length(v))
            %  Construct u for Householder reflector Pj+1.
            u = [zeros(initer,1); v(initer+1:end)];
            alpha = norm(u);
            if (alpha ~= 0)
                alpha = scalarsign(v(initer+1))*alpha;
                %  u = v(initer+1:end) +
                %        sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
                u(initer+1) = u(initer+1) + alpha;
                u = u / norm(u);
                U(:,initer+1) = u;
                
                %  Apply Pj+1 to v.
                %  v = v - 2*u*(u'*v);
                v(initer+2:end) = 0;
                v(initer+1) = -alpha;
            end
        end
        
        %  Apply Given's rotations to the newly formed v.
        for colJ = 1:initer-1
            tmpv = v(colJ);
            v(colJ)   = conj(J(1,colJ))*v(colJ) + conj(J(2,colJ))*v(colJ+1);
            v(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*v(colJ+1);
        end
        
        %  Compute Given's rotation Jm.
        if ~(initer==length(v))
            rho = norm(v(initer:initer+1));
            J(:,initer) = v(initer:initer+1)./rho;
            w(initer+1) = -J(2,initer).*w(initer);
            w(initer) = conj(J(1,initer)).*w(initer);
            v(initer) = rho;
            v(initer+1) = 0;
        end
        
        R(:,initer) = v(1:inner);
        
        normr = abs(w(initer+1));
        resvec((outiter-1)*inner+initer+1) = normr;
        normr_act = normr;
        
        if (normr <= tolb || stag >= maxstagsteps || moresteps)
            if evalxm == 0
                ytmp = R(1:initer,1:initer) \ w(1:initer);
                additive = U(:,initer)*(-2*ytmp(initer)*conj(U(initer,initer)));
                additive(initer) = additive(initer) + ytmp(initer);
                for k = initer-1 : -1 : 1
                    additive(k) = additive(k) + ytmp(k);
                    additive = additive - U(:,k)*(2*(U(:,k)'*additive));
                end
                if norm(additive) < eps*norm(x)
                    stag = stag + 1;
                else
                    stag = 0;
                end
                xm = x + additive;
                evalxm = 1;
            elseif evalxm == 1
                addvc = [-(R(1:initer-1,1:initer-1)\R(1:initer-1,initer))*...
                    (w(initer)/R(initer,initer)); w(initer)/R(initer,initer)];
                if norm(addvc) < eps*norm(xm)
                    stag = stag + 1;
                else
                    stag = 0;
                end
                additive = U(:,initer)*(-2*addvc(initer)*conj(U(initer,initer)));
                additive(initer) = additive(initer) + addvc(initer);
                for k = initer-1 : -1 : 1
                    additive(k) = additive(k) + addvc(k);
                    additive = additive - U(:,k)*(2*(U(:,k)'*additive));
                end
                xm = xm + additive;
            end
            r = b - iterapp('mtimes',afun,atype,afcnstr,xm,varargin{:});
            if norm(r) <= tol*n2b
                x = xm;
                flag = 0;
                iter = [outiter, initer];
                break
            end
            minv_r = r;
            if existM1
                
                %minv_r = P*minv_r;
                minv_r = r-A*(ZH*(AH\(ZH'*r)));
                minv_r = iterapp('mldivide',m1fun,m1type,m1fcnstr,minv_r,varargin{:});
                if ~all(isfinite(minv_r))
                    flag = 2;
                    break
                end
            end
            if existM2
                minv_r = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_r,varargin{:});
                if ~all(isfinite(minv_r))
                    flag = 2;
                    break
                end
            end
            
            normr_act = norm(minv_r);
            resvec((outiter-1)*inner+initer+1) = normr_act;
            
            if normr_act <= normrmin
                normrmin = normr_act;
                imin = outiter;
                jmin = initer;
                xmin = xm;
                minupdated = 1;
            end
            
            if normr_act <= tolb
                x = xm;
                flag = 0;
                iter = [outiter, initer];
                break
            else
                if stag >= maxstagsteps && moresteps == 0
                    stag = 0;
                end
                moresteps = moresteps + 1;
                if moresteps >= maxmsteps
                    if ~warned
                        warning(message('MATLAB:gmres:tooSmallTolerance'));
                    end
                    flag = 3;
                    iter = [outiter, initer];
                    break;
                end
            end
        end
        
        if normr_act <= normrmin
            normrmin = normr_act;
            imin = outiter;
            jmin = initer;
            minupdated = 1;
        end
        
        if stag >= maxstagsteps
            flag = 3;
            break;
        end
    end         % ends inner loop
    
    evalxm = 0;
    
    if flag ~= 0
        if minupdated
            idx = jmin;
        else
            idx = initer;
        end
        y = R(1:idx,1:idx) \ w(1:idx);
        additive = U(:,idx)*(-2*y(idx)*conj(U(idx,idx)));
        additive(idx) = additive(idx) + y(idx);
        for k = idx-1 : -1 : 1
            additive(k) = additive(k) + y(k);
            additive = additive - U(:,k)*(2*(U(:,k)'*additive));
        end
        x = x + additive;
        xmin = x;
        r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        minv_r = r;
        if existM1
            

            minv_r = r-A*(ZH*(AH\(ZH'*r)));
            minv_r = iterapp('mldivide',m1fun,m1type,m1fcnstr,minv_r,varargin{:});
            if ~all(isfinite(minv_r))
                flag = 2;
                break
            end
        end
        if existM2
            minv_r = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_r,varargin{:});
            if ~all(isfinite(minv_r))
                flag = 2;
                break
            end
        end
        normr_act = norm(minv_r);
        r = minv_r;
    end
    
    if normr_act <= normrmin
        xmin = x;
        normrmin = normr_act;
        imin = outiter;
        jmin = initer;
    end
    
    if flag == 3
        break;
    end
    if normr_act <= tolb
        flag = 0;
        iter = [outiter, initer];
        break;
    end
    minupdated = 0;
end         % ends outer loop

% returned solution is that with minimum residual
if flag == 0
    relres = normr_act / n2minv_b;
else
    x = xmin;
    iter = [imin jmin];
    relres = normr_act / n2minv_b;
end

% truncate the zeros from resvec
if flag <= 1 || flag == 3
    resvec = resvec(1:(outiter-1)*inner+initer+1);
    indices = resvec==0;
    resvec = resvec(~indices);
else
    if initer == 0
        resvec = resvec(1:(outiter-1)*inner+1);
    else
        resvec = resvec(1:(outiter-1)*inner+initer);
    end
end

% only display a message if the output flag is not used
% if nargout < 2
%     if restarted
%         itermsg(sprintf('gmres(%d)',restart),tol,maxit,[outiter initer],flag,iter,relres);
%     else
%         itermsg(sprintf('gmres'),tol,maxit,initer,flag,iter(2),relres);
%     end
% end

function sgn = scalarsign(d)
sgn = sign(d);
if (sgn == 0)
    sgn = 1;
end

function y = iterapp(op,afun,atype,afcnstr,x,varargin)
%ITERAPP   Apply matrix operator to vector and error gracefully.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X) applies matrix operator AFUN to vector
%   X. If ATYPE is 'matrix, then AFUN is a matrix and the OP is applied
%   directly. OP is either 'mtimes' or 'mldivide'.
%   ATYPE and AFCNSTR are used in case of error.
%   ITERAPP(OP,AFUN,ATYPE,AFCNSTR,X,P1,P2,...) allows extra arguments to
%   AFUN(X,P1,P2,...) although this usage is now discouraged in favor of
%   using anonymous functions.
%   AFUN(X,P1,P2,...,PN,TFLAG) should accept a TFLAG as its final input if
%   the calling function is BICG, LSQR or QMR. TFLAG is either 'transp' or
%   'notransp' depending on whether A' OP X or A OP X is required.
%   ITERAPP is designed for use by iterative methods like PCG which
%   require matrix operators AFUN representing matrices A to operate on
%   vectors X and return A*X and may also have operators MFUN representing
%   preconditioning matrices M operate on vectors X and return M\X.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.7.4.8 $ $Date: 2011/05/17 02:33:16 $

if strcmp(atype,'matrix')
    switch lower(op)
        case 'mtimes'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' * x;
            else
                y = afun * x;
            end
        case 'mldivide'
            if (nargin >= 6) && isequal(varargin{end},'transp')
                y = afun' \ x;
            else
                y = afun \ x;
            end
        otherwise
            error(message('MATLAB:iterapp:InvalidOp'))
    end
else
    try
        if (nargin >= 6) && isequal(varargin{end},'notransp')
            % A request for A*x coming from BICG, LSQR and QMR
            try
                % New syntax: we now request afun(x,P1,P2,...,PN,'notransp')
                y = afun(x,varargin{:});
            catch
                % Old syntax: we used to accept afun(x,P1,P2,...,PN)
                y = afun(x,varargin{1:end-1});
            end
        else
            % A request for A*x
            % coming from BICGSTAB, CGS, GMRES, MINRES, PCG or SYMMLQ
            % with the call always afun(P1,P2,...,PN)
            % or a request for A'*x coming from
            % BICG, LSQR and QMR in the afun(x,P1,P2,...,PN,'transp') case
            y = afun(x,varargin{:});
        end
    catch ME
        error(message('MATLAB:iterapp:InvalidInput', atype,afcnstr, ME.message));
    end

    if ~iscolumn(y)
        error(message('MATLAB:iterapp:MustReturnColumn', atype, afcnstr));
    end
end


function [atype,afun,afcnstr] = iterchk(A)
%ITERCHK  Checks arguments to iterative methods.
%   [ATYPE,AFUN,AFCNSTR] = ITERCHK(A) returns the following:
%   ATYPE is either 'matrix', 'function', 'expression' or 'inline object'.
%   AFUN is the function name or inline object.
%   AFUN is '' if ATYPE is 'matrix'.
%   AFCNSTR is the function name if ATYPE is 'function'.
%   AFCNSTR is the formula of the function if ATYPE is 'expression' or
%   'inline object'.  AFCNSTR is '' if ATYPE is 'matrix'.
%
%   See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 1.8.4.3 $ $Date: 2010/08/23 23:12:57 $


[afun,afunmsg] = fcnchk(A);
if isempty(afunmsg)
   if isa(afun,'inline')      
      if isa(A,'inline')
         atype = 'inline object';
      else
         atype = 'expression';
      end
      afcnstr = formula(afun);
   else % both function_handles @fun and function names 'fun'
      atype = 'function';
      if isa(A,'function_handle')
          afcnstr = func2str(A);
      else
          afcnstr = A;
      end
   end
elseif isa(A,'float')
   afun = A;
   atype = 'matrix';
   afcnstr = '';
else
   error(message('MATLAB:iterchk:InvalidInput'));
end


























