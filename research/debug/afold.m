function [R,S] = afold(A,As,isM)

%% Index
N = size(A,1);
isF = isM;   % fine node
isC = ~isF;  % coarse node
Nf = sum(isF);

%% Approximated factorization
Dfsqrt = sqrt(full(diag(A(isF,isF))));  % diag of Af
Dinv = spdiags(1./Dfsqrt,0,Nf,Nf);
Afc = Dinv*A(isF,isC);

%% Append Af to R
R = sparse(Nf,N);
R(1:Nf,1:Nf) = spdiags(Dfsqrt,0,Nf,Nf);
R(1:Nf,Nf+1:N) = Afc;

%% Approximated Schur complement
% S = A(isC,isC) - A(isC,isF)*inv(A(isF,isF))*A(isF,isC);
[Pro,Res] = interpolationAMGt(As,isC);
% [Pro,Res] = interpolationAMGa(As,isC);
% [Pro,Res] = interpolation4achol(A,isC);
% [Pro,Res] = interpolationAMGs(As,isC);
S = Res*A*Pro;