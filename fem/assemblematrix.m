function [A,M,area] = assemblematrix(node,elem,lumpflag,K)
%% ASSEMBLEMATRIX matrix for diffusion and reaction
%
% [A,M] = ASSEMBLEMATRIX(node,elem) return the stiffness matrix and the
% mass matrix.
%
% [A,M] = ASSEMBLEMATRIX(node,elem,1) return the stiffness matrix and the
% lumped mass matrix. Note that in the output M is a vector not a matrix. A
% sparse diagonal matrix using M as diaongal can be obtained by
% spdiags(M,0,N,N);
%
% A = ASSEMBLEMATRIX(node,elem,[],K) returns stiffness matrix with
% piecewise constant coefficient K.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Parameters
N = size(node,1);
A = sparse(N,N);
M = sparse(N,N);
if ~exist('K','var'), K = []; end
if ~exist('lumpflag','var'), lumpflag = 0; end
if (nargout > 1)
    if ~exist('lumpflag','var') || isempty(lumpflag)
        lumpflag = 0;
    end    
end

%% 3-D case
if (size(node,2) == 3) && (size(elem,2) == 4) % 3-D 
    [A,M,area] = assemblematrix3(node,elem,lumpflag);
    return
end

%% Compute vedge, edge as a vector, and area of each element
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));

%% Assemble stiffness matrix
for i = 1:3
    for j = 1:3
        Aij = (ve(:,1,i).*ve(:,1,j)+ve(:,2,i).*ve(:,2,j))./(4*area);
        if ~isempty(K), Aij = K.*Aij; end
        A = A + sparse(elem(:,i),elem(:,j),Aij,N,N);
        if ~lumpflag 
           Mij = area*((i==j)+1)/12;
           M = M + sparse(elem(:,i),elem(:,j),Mij,N,N);
        end
    end
end

%% Assemble the mass matrix by the mass lumping
if lumpflag
    M = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
end