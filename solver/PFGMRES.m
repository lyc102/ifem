function [u, iter, error] = PFGMRES(A, f, u, maxit, restart, tol, M, print_level)
% Preconditioned Flexible General Minimal Residual Method (Right Preconditioner)
%
% BRIEF: 
%   Solving Au=f by general minimal residual using flexible preconditioner
%   M (right preconditioning)
%
%
% INPUT: 
%   A:              Matrix (can be a function handle)
%   f:              Right hand side
%   u:              Initial guess
%   maxit:          Maximial number of itrations allowed
%   restart:        Restart number
%   tol:            Tolerance
%   M:              Precodnitoner (can be a function handle)
%   print_level:    How much information to print (=0: no output; >0 output)
%
% OUTPUT:
%   u:              Solution
%   iter:              Number of iterations
%   error:          History of l^2 norm of residual
%
% USAGE:
%   [u, iter, error] = PFGMRES(A, f, u, maxit, restart, tol, [],
%       print_level): PFGMRES without preconditoner
%   [u, iter, error] = PFGMRES(A, f, u, maxit, restart, tol, M,
%       print_level): PFGMRES using preconditoner M
%
% COPYRIGHT:
% 
%   X.Hu 02/03/2011, Penn State University

% TODO:
%   1. Check false convergence

%-------------------
% Preparation 
%-------------------
% size of the problem
N = size(f,1);

% restart number
if restart ~= 0
    restart = min(restart, N);
    restart = min(restart,maxit);
else
    restart = min(maxit,N);
end

% initalize memory
r = zeros(N,1);
V = zeros(N, restart+1);
Z = zeros(N, restart);
H = zeros(restart+1, restart);
b = zeros(restart+1,1);
R = zeros(restart,restart);
c = zeros(restart,1);
s = zeros(restart,1);
y = zeros(restart,1);
error = zeros(maxit+1,1);

%local variables
iter = 1;
%norm_r = 0.0;

% computer the residual
if isa(A, 'double')
    r = f - A*u;
elseif isa(A, 'function_handle')
    r = f - A(u);
else
    error('A is neither a matrix or a function handle!!!');
end % end if

% store the residual
normr = norm(r);
error(1) = normr;

% output if needed
if (print_level > 0)
    fprintf('#It|  ||r||/||r0||  |  ||r||   |    Conv. Factor |\n');
    fprintf(' %d |  %e  |  %e  | %f |\n', 0, 1.0, error(1), 0.0);
end

%-------------------
% Main loop 
%-------------------
while (iter < maxit)
    
    % reset converge  
    converge = 0;

    % first orthonormal basis v1
    norm_r = norm(r);
    V(:,1) = r/norm_r;
    
    % form right hand side b for the hessenberg system
    b(1) = norm_r; 
    
    % loop for restart
    for i = 1:restart
        
        % Preconditioning: z = M\v
        if isempty(M)
            Z(:,i) = V(:,i);
        elseif isa(M, 'double')
            Z(:,i) = M\V(:,i);
        elseif isa(M, 'function_handle')
            Z(:,i) = M(V(:,i));
        else
            error('Preconditoner M is invalid!!!');
        end % end if
        
        % w = Az
        if isa(A, 'double')
            V(:,i+1) = A*Z(:,i);
        elseif isa(A, 'function_handle')
            V(:,i+1) = A(Z(:,i));
        else
            error('A is neither a matrix or a function handle!!!');
        end % end if
        
        %--------------------------------------
        % modified Gram-Schmidt 
        %--------------------------------------
        for k = 1:i
            
          H(k,i) = V(:,k)'*V(:,i+1);
          V(:,i+1) = V(:,i+1) - H(k,i)*V(:,k);
            
        end % end for k
        
        % new orthonormal basis 
        H(i+1,i) = norm(V(:,i+1));
        V(:,i+1) = V(:,i+1)/H(i+1,i); % becareful small H(i+1,i)
        
        %--------------------------------------
        % Use Givens transformation to get upper triangular system R 
        %--------------------------------------
        R(1,i) = H(1,i);
        
        % apply the previous Givens transformations
        if (i~=1)
            
            for k = 2:i
               
                temp = c(k-1)*R(k-1,i) + s(k-1)*H(k,i);
                R(k,i) = -s(k-1)*R(k-1,i) + c(k-1)*H(k,i);
                R(k-1,i) = temp;
                
            end % end for k
            
        end % end if (i~=1)
        
        % new Givens transformation
        delta = sqrt(R(i,i)^2 + H(i+1,i)^2);
        c(i) = R(i,i)/delta;
        s(i) = H(i+1,i)/delta;
        
        R(i,i) = c(i)*R(i,i) + s(i)*H(i+1,i);
        
        % apply Givens transformation to Right hand side b
        b(i+1) = -s(i)*b(i);
        b(i) = c(i)*b(i);
        
        % count iterations
        iter = iter + 1;
        
        % check convergence b(i+1) = || f-Au_k ||  
        error(iter) = abs(b(i+1));
        
        % output
        if (print_level > 0)
            fprintf(' %d |  %e  |  %e  | %f |\n', iter-1, error(iter)/error(1), error(iter), error(iter)/error(iter-1));
        end
        
        if ((error(iter)/error(1)) < tol)
            converge = 1;
            break;
        end
        
        %--------------------------------------
        
    end % end for i
    
    % solve the upper trangular matrix
    y(1:i) = R(1:i, 1:i)\b(1:i);
    
    % solution
    u = u + Z(:,1:i)*y(1:i);
    
    % check convergence
%     if (converge)
%        break;
%     end
    
    % update residual and restart
    if isa(A, 'double')
        r = f - A*u;
    elseif isa(A, 'function_handle')
        r = f - A(u);
    else
        error('A is neither a matrix or a function handle!!!');
    end % end if
    
    % check convergence
    if (converge && ((norm(r)/normr) < tol))
        break;
    end
    
end % end while iter

end

