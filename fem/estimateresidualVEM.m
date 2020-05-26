function [eta, etaTotal, elemResTotal] = estimateresidualVEM(node,elem,u,pde)
%% ESTIMATERESIDUALVEM residual type error estimator for 2D virtual element.
%
% eta = estimateresidualVEM(node,elem,u,pde) computes the residual type a
% posterior error estimator for the Poisson equation in the form
%
% $$\eta _T^2 = \sum _{E\in T}[Pi_E \nabla u\cdot v_E^{\bot}]^2 
%               + area(T)^2\sum _i w_i f_i^2,$$
%
% which is an approximation of 
%
% $$\eta _T^2 =\sum _{E\in T}\int _E h_E[\nabla u\cdot n_E]^2 ds 
%               + \int _T h_E^2 f^2 dx.$$
%
% Here [Pi_E \nabla u\cdot n_E] denotes the jump of the projected flux accross the
% interior edge E and n_E is the normal vector of E. In 2-D, n_E is the
% right 90 degree rotation of the edge vector of E assuming E is oriented 
% counter-clockwisely. 
% The integral $\int_E h_E[Pi_E\nabla u\cdot n_E]^2 ds$ is simplified to 
% $[Pi_E \nabla u\cdot v_E^{\bot}]^2$.
%
% For general diffusion equation div(K grad u) = f, the flux is K grad u.
% The jump will be [K Pi_E\nabla u\cdot n_E].
%
% PDE is a structure array that records the information of Neumann
% boundary condition and diffusion coefficent.
%
%  pde.f     : right hand side
%  pde.d     : diffusion coeffcients
%
% TODO: Neumann boundary
% 
% Example
%
%   KelloggVEM
%
% See also estimateresidual, estimaterecoveryVEM
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Compute the flux from the finite element function u
N = size(node,1);
NT = size(elem,1);
[Du, area, center] = graduVEM(node,elem,u);


if exist('pde','var') && isfield(pde,'d') && ~isempty(pde.d)
    if isreal(pde.d)
        K = pde.d;                  % d is an array
    else                            % d is a function
        K = pde.d(center);              
    end
else
    K = ones(size(elem,1),1);
end
Du = [K.*Du(:,1),K.*Du(:,2)];   % flux

%% If the PDE is the reaction diffusion, modify the right hand side
if ~isfield(pde,'f') || (isreal(pde.f) && (~any(pde.f)))
    pde.f = [];
end

if ~isfield(pde,'c'), pde.c = []; end

if ~isempty(pde.c) % reaction constant C (piecewise linear preferred)
    if isnumeric(pde.c)
        C = pde.c;    % C is an array
    end
    if  ~isnumeric(pde.c) % C is a function
        C = pde.c(node); % make C a piecewise linear
    end
    if  isreal(pde.f) % f is the an array
        switch length(pde.f)
            case N % f is piecewise linear
                pde.f = pde.f - C.*u;
            case NT % f is piecewise constant
                % change element-wise value to node-wise value by averaging
                % on vertex patch
                v2t = cell2node(elem)';
                patchArea = v2t*area;
                f2V = v2t*pde.f./patchArea;
                pde.f = f2V - C.*u;      
        end
    elseif ~isreal(pde.f) % f is a function
        pde.f = pde.f(node)- C.*u;
    elseif isempty(pde.f) % f = 0
        pde.f = -C.*u;
    end
end
%% Jump of normal flux
T = auxstructurepoly(elem);
neighbor = T.neighbor; 
elemVertNum = cellfun('length',elem);% the number of vertices per element
minNv = min(elemVertNum);
maxNv = max(elemVertNum);
clear T
edgeJump = zeros(NT,1);
for nV = minNv:maxNv
    isNv = (elemVertNum == nV);
    NelemNv = sum(isNv);
    if NelemNv==0; continue; end
    
    elemNv = cell2mat(elem(isNv));
    neighborNv = cell2mat(neighbor(isNv));
    DuNv = Du(isNv,:);
    KnV = K(isNv,:);
    % edge vector
    x1 = reshape(node(elemNv,1),[NelemNv,nV]);
    y1 = reshape(node(elemNv,2),[NelemNv,nV]);
    x2 = circshift(x1,[0,-1]);
    y2 = circshift(y1,[0,-1]);
    normVecx = y2 - y1; % scaled normal vector
    normVecy = x1 - x2; % \nu_e = n_e |e|
    ne = permute(cat(3, normVecx, normVecy), [1 3 2]);
    
    for i = 1:nV
        edgeJump(isNv,:) = edgeJump(isNv,:) + ...
            dot((DuNv-Du(neighborNv(:,i),:)),ne(:,:,i),2).^2 ...
            ./(0.5*KnV+0.5*K(neighborNv(:,i)));
    end
end

%% Elementwise residual
elemResidual = zeros(size(elem,1),1);

if isfield(pde,'f') && ~isempty(pde.f) && ~isnumeric(pde.f)
    pde.f = pde.f(node);
end
if isnumeric(pde.f) && (length(pde.f) == size(elem,1))
    % piecewise constant
    elemResidual = pde.f.^2./K;
elseif isnumeric(pde.f) && (length(pde.f) == size(node,1))
    % nodal
    f2elem = pde.f([elem{:}]');
    f2elem = array2cell(f2elem, elemVertNum);
    elemResidual = cellfun(@mean, f2elem).^2./K;
end

elemResidual = elemResidual.*area.^2;
elemResidual(isnan(elemResidual)) = 0; % singular values are excluded
%% Residual type error estimator
eta = (elemResidual + edgeJump).^(1/2);
etaTotal = norm(eta);
elemResTotal = sqrt(sum(elemResidual));