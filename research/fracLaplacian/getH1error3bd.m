function err = getH1error3bd(node,Neumann,pde,u,A,quadorder)

if ~exist('quadorder','var'), quadorder = 4; end

v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
v14 = node(Neumann(:,4),:)-node(Neumann(:,1),:);
area = sqrt(abs(sum(mycross(v12,v14,2).^2,2)));

NT = size(Neumann,1);
err = zeros(NT,1);
[GaussPt,weightgN] = quadptsquad(quadorder);
lambdax1 = GaussPt(:,1);
lambdax2 = GaussPt(:,2);
nQuadgN = size(GaussPt,1);
for pp = 1:nQuadgN
    phi = kron([lambdax2(pp) 1-lambdax2(pp)],[lambdax1(pp); 1-lambdax1(pp)]);
    phi = phi([1 2 4 3]); % switch the index of 3 and 4
    % quadrature points in the x-y coordinate
    ppxyz = node(Neumann(:,1),:)*phi(1) + node(Neumann(:,2),:)*phi(2)...
          + node(Neumann(:,3),:)*phi(3) + node(Neumann(:,4),:)*phi(4);
    gNp = pde.g_N(ppxyz);    
    uhp = u(Neumann(:,1))*phi(1) + u(Neumann(:,2))*phi(2) ...
        + u(Neumann(:,3))*phi(3) + u(Neumann(:,4))*phi(4);
    up = pde.exactu(ppxyz);
    err = err + weightgN(pp)*gNp.*(up - 2*uhp);
end
err = sum(err.*area) + u'*A*u;
err = sqrt(abs(err));
