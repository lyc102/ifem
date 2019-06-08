function [Pro,FreeNodec] = transferP1red(elemc,elemf,FreeNode)
%% TRANSFERP2RED
%
%

if ~exist('FreeNode','var'), FreeNode = []; end 

%% Data structure
N = max(elemc(:));
if exist('elemf','var')
    [elem,HB] = uniformcoarsenred(elemf);
else
    [elem2dof,edge] = dofP2(elemc);
    NE = size(edge,1);
    HB(:,[1 2 3]) = [(N+1:N+NE)', edge(:,1:2)]; 
end

%% Prolongation matrix
coarseNode = (1:N)';
