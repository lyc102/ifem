function [QDu,area,center,QDlambda,stbElem] = graduVEM(node,elem,u)
%% GRADUVEM the projected gradient of a linear virtual element function.
%
% QDu = GRADUVEM(node,elem,u) compute the L2 projection of the gradient of 
% a virtual element function u on a polygonal mesh representing by (node,elem).
% 
% [QDu,area,center,QDlambda,stbElem] = GRADUVEM(node,elem,u) also outputs areas, centroids, 
% elementwise stabilization, and the L2 projection of Dlambda which is the gradient of P1 
% conforming VEM basis. QDu{t}(i,1) is the x-component of the i-th vertex's basis in t-th element.
% stbElem(t,:) is the \ell^2 norm of (I-P)u_h on each element.
% 
% Remark: this routine is mostly following Long's vectorization in PoissonVEM,
% however is ONLY implemented for Matlab version > R2016b
% see: https://scaomath.github.io/blog/MATLAB-update-on-array-compatibility/
%
% See also gradu, gradbasis, PoissonVEM
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.
%%
if ~iscell(elem); elem = num2cell(elem,2); end

%%
NT = size(elem,1);
QDu = zeros(NT,2);
area = zeros(NT,1);
center = zeros(NT,2);
QDlambda = cell(NT,1);
stbElem = zeros(NT,1); % || (I- P) u ||_{l^2,K}
%%
elemVertexNumber = cellfun('length',elem);% the number of vertices per element
minNv = min(elemVertexNumber);
maxNv = max(elemVertexNumber);


%% 
for nV = minNv:maxNv
    isNv = (elemVertexNumber == nV); % index of elements with Nv vertices
    if ~any(isNv); continue; end
    elemNv = cell2mat(elem(isNv));
    NelemNv = sum(isNv);
    x1 = reshape(node(elemNv,1),[NelemNv,nV]);
    y1 = reshape(node(elemNv,2),[NelemNv,nV]);
    x2 = circshift(x1,[0,-1]);
    y2 = circshift(y1,[0,-1]);
    bdIntegral = x1.*y2 - y1.*x2;
    areaNv = sum(bdIntegral,2)/2; 
    xc = sum((x1+x2).*bdIntegral,2)./abs(areaNv)/6;
    yc = sum((y1+y2).*bdIntegral,2)./abs(areaNv)/6;
    center(isNv,:) = [xc, yc];
    normVecx = y2 - y1; % normal vector is a rotation of edge vectors
    normVecy = x1 - x2;
    Bx = (normVecx + circshift(normVecx,[0,1]))./(2*areaNv); % average of normal vectors
    By = (normVecy + circshift(normVecy,[0,1]))./(2*areaNv); % in adjacent edges
    QDlambdaNv = reshape([Bx, By]',[nV,2,NelemNv]);
    QDlambda(isNv) = squeeze(num2cell(QDlambdaNv, [1, 2]));
    uh2elemNv = u(elemNv);
    if NelemNv == 1; uh2elemNv = uh2elemNv'; end % in this case: size(uh2elemNv,2) == 1
    QDudx = sum(uh2elemNv.*Bx, 2);
    QDudy = sum(uh2elemNv.*By, 2);
    QDlambda(isNv) = squeeze(num2cell(QDlambdaNv, [1, 2]));
    %%
    area(isNv,:) = abs(areaNv);
    QDu(isNv,:) = [QDudx, QDudy];
    
    %% compute stablization
    if nargout > 4
        h = sign(areaNv).*sqrt(abs(areaNv)); % h = sqrt(area) not the diameter
        cx = mean(reshape(node(elemNv(:),1),[NelemNv,nV]),2);
        cy = mean(reshape(node(elemNv(:),2),[NelemNv,nV]),2);
        Dx = (x1 - cx)./h; Dy = (y1 - cy)./h; %  m = (x - cx)/h
        IminusP = zeros(NelemNv,nV,nV);
        for i = 1:nV
            for j = 1:nV
                IminusP(:,i,j) = - 1/nV - Dx(:,i).*Bx(:,j).*abs(h) - Dy(:,i).*By(:,j).*abs(h);
            end
            IminusP(:,i,i) = ones(NelemNv,1) + IminusP(:,i,i);
        end
        if NelemNv == 1
            IminusP = squeeze(IminusP);
            utIminusPu = uh2elemNv*IminusP*uh2elemNv';
        else
            IminusPu = sum(bsxfun(@times, IminusP, uh2elemNv), 2);
            utIminusPu = sum(uh2elemNv.*squeeze(IminusPu),2);
        end
        
        stbElem(isNv) = sqrt(abs(utIminusPu));
    end
end