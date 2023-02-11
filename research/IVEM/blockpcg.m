function [x,info] = blockpcg(A,b,x0,tol,maxIt,D,L,U,P,Q)


k = 1;
% if printlevel>=1
%     fprintf('Conjugate Gradient Method using HX preconditioner \n');
% end
x = x0;
r = b - A*x;
nb = norm(b);
err = zeros(maxIt,2);
err(1,:) = norm(r)/nb;
blkId = length(D);
NumDof = size(b,1);

while (max(err(k,:)) > tol) && (k <= maxIt)
    % compute Br by HX preconditioner
    r1 = D.^(-1).*r(1:blkId);
    r2 = P*r(blkId+1:NumDof);
    r2 = L\r2;
    r2 = U\r2;
    r2 = Q*r2;
    Br = [r1;r2];
    
    % update tau, beta, and p
    rho = r'*Br;  % r'*Br = e'*ABA*e approximates e'*A*e
    if k==1
        p = Br;
    else
        beta = rho/rho_old;
        p = Br + beta*p;
    end
    % update alpha, x, and r
    Ap = A*p;
    alpha = rho/(Ap'*p);
    r = r - alpha*Ap;
    x = x + alpha*p;
    rho_old = rho;
    % compute err for the stopping criterion
    k = k + 1;
    err(k,1) = sqrt(abs(rho/(x'*b))); % approximate relative error in energy norm
    err(k,2) = norm(r)/nb; % relative error of the residual in L2-norm
%     if printlevel >= 2
%         fprintf('#dof: %8.0u, BCG iter: %2.0u, err = %12.8g\n',...
%             Ndof, k, max(err(k,:)));
%     end
end
err = err(1:k,:);
itStep = k-1;
if k > maxIt || (max(err(end,:))>tol)
    flag = 1;
else
    flag = 0;
end

info = struct('itStep',itStep,'error',err,'flag',flag,'stopErr',max(err(end,:)));   