function M = getmassmatrixP2(elem2dof,area,type)
%% GETMASSMATRIXP2  get the mass matrix of P2 elements
%
% Type:
% - HB:
% - NB:
% - NBB:
%
% Added by Lin Zhong. Improved by Long Chen.
%
%% Data Structure
NT = length(area);
Ndof = double(max(elem2dof(:)));

%% Mass matrices for NBB bases
if strcmp(type,'NBB')
    elem2dof(:,7) = uint32(Ndof+(1:NT)'); % add bubble index    
    Ndof = double(max(elem2dof(:)));
    elemM(1:NT,7)   = 9/20;
    elemM(1:NT,1:3) = 1/20;
    elemM(1:NT,4:6) = 2/15;
    elemM = elemM.*repmat(area,1,7);
    diagM = accumarray(elem2dof(:), elemM(:), [Ndof 1]);
    M = spdiags(double(diagM),0,Ndof,Ndof);
    return;
end

%% Mass matrices for HB or NB bases
[lambda, w] = quadpts(4);
nQuad = size(lambda,1);
if strcmp(type,'HB')
    phi(:,1) = lambda(:,1);
    phi(:,2) = lambda(:,2);
    phi(:,3) = lambda(:,3);
elseif strcmp(type,'NB')
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
end
phi(:,4) = 4*lambda(:,2).*lambda(:,3);
phi(:,5) = 4*lambda(:,3).*lambda(:,1);
phi(:,6) = 4*lambda(:,1).*lambda(:,2);

Nbases = 6; NMhalf = 21;
ii = zeros(NMhalf*NT,1); 
jj = zeros(NMhalf*NT,1); 
sM = zeros(NMhalf*NT,1);
index = 0;
for i = 1:Nbases
    for j = i:Nbases
        Mij = 0;
        for p = 1:nQuad
            Mij = Mij + w(p)*phi(p,i).*phi(p,j);
        end
        Mij = Mij.*area;
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));
        sM(index+1:index+NT) = Mij;
        index = index + NT;
    end
end
clear Mij
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
M = M + MU + MU';