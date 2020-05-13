function [QDu,area,QDlambda] = graduVEM(node,elem,u)
%% GRADUVEM the projected gradient of a linear virtual element function.
%
% QDu = GRADUVEM(node,elem,u) compute the L2 projection of the gradient of 
% a virtual element function u on a polygonal mesh representing by (node,elem).
% 
% [QDu,area,QDlambda] = GRADUVEM(node,elem,u) also outputs area and L2 projection 
% of Dlambda which is the gradient of P1 conforming VEM basis. QDu{t}(i,1) is
% the x-component of the i-th vertex's basis in t-th element.
%
% See also gradu, gradbasis
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.
%%
if ~iscell(elem); elem = num2cell(elem,2); end

%%
NT = size(elem,1);
QDu = zeros(NT,2);
area = zeros(NT,1);
QDlambda = cell(NT,1);
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
    normVecx = y2 - y1; % normal vector is a rotation of edge vector
    normVecy = x1 - x2;
    Bx = (normVecx + circshift(normVecx,[0,1]))./(2*areaNv); % average of normal vectors
    By = (normVecy + circshift(normVecy,[0,1]))./(2*areaNv); % in adjacent edges
    QDlambdaNv = reshape([Bx, By]',[nV,2,NelemNv]);
    QDlambda(isNv) = squeeze(num2cell(QDlambdaNv, [1, 2]));
    QDudx = sum(u(elemNv).*Bx, 2);
    QDudy = sum(u(elemNv).*By, 2);
    QDlambda(isNv) = squeeze(num2cell(QDlambdaNv, [1, 2]));
    %%
    area(isNv,:) = areaNv;
    QDu(isNv,:) = [QDudx, QDudy];
end