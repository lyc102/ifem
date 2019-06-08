function err = getH1errorbd(node,Neumann,pde,u,A,quadorder)

if ~exist('quadorder','var'), quadorder = 4; end

el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
[lambdagN,weightgN] = quadpts1(quadorder);
phigN = lambdagN;                 % linear bases
nQuadgN = size(lambdagN,1);
err = zeros(size(Neumann,1),1);
for pp = 1:nQuadgN
    % quadrature points in the x-y coordinate
    ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
         + lambdagN(pp,2)*node(Neumann(:,2),:);
    uhp = u(Neumann(:,1))*phigN(pp,1) + u(Neumann(:,2))*phigN(pp,2);
    up = pde.exactu(ppxy);
    gNp = pde.g_N(ppxy);
    err = err + weightgN(pp)*gNp.*(up - 2*uhp);
end
err = sum(err.*el) + u'*A*u;
err = sqrt(abs(err));