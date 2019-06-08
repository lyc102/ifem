function [i,j,s,Ac] = af(A,As,isM)
%%
% i,j: index system of A

%% Index
% N = size(A,1);
F = find(isM);  % fine node
isC = ~isM;
C = find(isC); % coarse node
Nf = length(F);

%% Approximated factorization
Dfsqrt = sqrt(full(diag(A(F,F))));  % diag of Af
Dinv = spdiags(1./Dfsqrt,0,Nf,Nf);
Acf = A(C,F)*Dinv;
% Afc = Dinv*A(F,isC);

%% Store L using coordinate format
[iC,jF,s] = find(Acf);  % i,j: index of C and F 
i = C(iC); 
j = F(jF);
i(end+1:end+Nf) = F;
j(end+1:end+Nf) = F;
s(end+1:end+Nf) = Dfsqrt;

%% Approximated Schur complement
if isempty(C)
    Ac = [];
    return
end
% Ac = A(isC,isC) - Afc'*Afc;
[Pro,Res] = interpolationAMGt(As,isC);
% [Pro,Res] = interpolationAMGa(As,isC);
% [Pro,Res] = interpolation4achol(A,isC);
% [Pro,Res] = interpolationAMGs(As,isC);
Ac = Res*A*Pro;