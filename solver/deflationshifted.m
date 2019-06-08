function [u,res,iter] = deflationshifted(elem,A,b,M,isfreeNode,tol)
global ZH AH;
[HB, NL, level]      = HBstructure(elem,10);
% isfreeNode           = false(size(node,1),1);
% isfreeNode(freeNode) = 1;
[Pro]                = transferoperator(HB,NL,isfreeNode);
ZH                   = Pro{level-1};
AH                   = ZH'*A*ZH;
%invAH = inv(AH);
%P = (speye(Ni)-AD(freeNode,freeNode)*ZH*invAH*ZH');
u1 = ZH*(AH\(ZH'*b));
[u,~,res,iter]  = Pgmres_deflation(A,b,[],tol,size(b,1),M);
iter = max(iter);
u2 = u - ZH*(AH\(ZH'*(A*u)));
u0 = u1+u2;
u  = u0;
end
















