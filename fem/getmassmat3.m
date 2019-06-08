function M = getmassmat3(node,elem2dof,volume,type,K)
%% GETMASSMAT Get mass matrix of the finite element space
%
% M = GETMASSMAT(node,elem2dof,volume,type,K) get mass matrix of the finite element
% space specified by elemType.  
%
% The type can be: 
% - P1: full mass matrix for P1 element
% - lump: lumped mass matrix for P1 element

N = size(node,1);
NT = length(volume);
Ndof = double(max(elem2dof(:)));

%% Coefficients and default type
if ~exist('type','var'), type = 'P1'; end
if ~exist('K','var'), K = []; end

%% Assembling
n = size(elem2dof,2);
switch n
    case 4  % P1 element
    if strcmp(type,'lump')
        %% Assemble the mass matrix by the mass lumping
        M = accumarray([elem2dof(:,1);elem2dof(:,2);elem2dof(:,3);elem2dof(:,4)],...
                       [volume;volume;volume;volume]/4,[N,1]);    
        if exist('K','var') && ~isempty(K) && ~isnumeric(K) % K is a function
            M = K(node).*M;
        elseif exist('K','var') && ~isempty(K) && isnumeric(K) && size(K,1) == N 
            M = K.*M;
        end 
        M = spdiags(M,0,N,N);                   
    else
        %% Assemble the full mass matrix
        if exist('K','var') && ~isempty(K) && ~isnumeric(K) % K is a function
            center = (node(elem2dof(:,1),:) + node(elem2dof(:,2),:) + ...
                      node(elem2dof(:,3),:) + node(elem2dof(:,4),:))/4;
            volume = K(center).*volume;
        elseif exist('K','var') && ~isempty(K) && isnumeric(K) && size(K,1) == NT 
            volume = K.*volume;
        end 
        M = sparse(N,N);
        for i = 1:4
            for j = i:4
                ii = double(elem2dof(:,i));
                jj = double(elem2dof(:,j));
                if (j==i)
                    M = M + sparse(ii,jj,volume/10,N,N);
                else
                    M = M + sparse([ii;jj],[jj;ii],[volume/20; volume/20],N,N);                               
                end                    
            end
        end        
    end
    case 10 % P2 element        
end