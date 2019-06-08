function [A,M,volume] = assemblematrix3(node,elem,lumpflag)
%% ASSEMBLEMATRIX3 matrix for diffusion and reaction
%
% [A,M] = ASSEMBLEMATRIX3(node,elem) return the stiffness matrix and the
% mass matrix.
%
% [A,M] = ASSEMBLEMATRIX3(node,elem,1) return the stiffness matrix and the
% lumped mass matrix. Note that in the output M is a vector not a matrix. A
% sparse diagonal matrix using M as diaongal can be obtained by
% spdiags(M,0,N,N);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = size(node,1);
A = sparse(N,N);
M = sparse(N,N);

%% Compute geometric quantities and gradient of local basis
[Dphi,volume] = gradbasis3(node,elem);

%% Assemble stiffness matrix
for i = 1:4
    for j = i:4
        Aij = dot(Dphi(:,:,i),Dphi(:,:,j),2).*volume;
        ii = double(elem(:,i));
        jj = double(elem(:,j));
        if (j==i)
            A = A + sparse(ii,jj,Aij,N,N);
        else
            A = A + sparse([ii;jj],[jj;ii],[Aij; Aij],N,N);        
        end        
        if (nargout > 1) && (~exist('lumpflag','var') || lumpflag == 0 )
            if (j==i)
                M = M + sparse(ii,jj,volume/10,N,N);
            else
                M = M + sparse([ii;jj],[jj;ii],[volume/20; volume/20],N,N);                               
            end                    
        end
    end
end

%% Assemble mass matrix by mass lumping
if exist('lumpflag','var') && lumpflag==1 && nargout > 1
    M = accumarray(elem(:),[volume;volume;volume;volume]/4,[N,1]);
end