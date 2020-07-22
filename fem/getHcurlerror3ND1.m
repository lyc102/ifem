function err = getHcurlerror3ND1(node,elem,curlE,Eh,markedElem)
%% GETHCURLERROR3ND1 Hcurl norm of approximation error for the linear Nedelect element in 3-D.
%
% err = GETHCURLERROR3ND1(node,elem,curlE,Eh,markedElem);
%
% NOTE: it is identical to getHcurlerror3ND since the added basis has no
% contribution to the curl part. 
%
% Example
% 
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     maxIt = 4;
%     HcurlErr = zeros(maxIt,1);
%     N = zeros(maxIt,1);
%     for k = 1:maxIt
%         [node,elem] = uniformbisect3(node,elem);
%         [elem2dof,edge] = dof3edge(elem);
%         pde = Maxwelldata2;
%         uI = edgeinterpolate1(pde.exactu,node,edge);
%         HcurlErr(k) = getHcurlerror3ND1(node,elem,pde.curlu,uI);
%         N(k) = length(uI);
%     end
%     r = showrate(N,HcurlErr,1,'b-+');
%     legend('||u-u_I||_{curl}',['N^{' num2str(r) '}'],'LOCATION','Best');
%
% See also getHcurlerror3ND1, getHcurlerror3ND2, getL2error3ND
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if nargin <=4
    err = getHcurlerror3ND(node,elem,curlE,Eh);
else
    err = getHcurlerror3ND(node,elem,curlE,Eh,markedElem);
end
%% TODO write more M-lint