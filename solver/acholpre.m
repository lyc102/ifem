function x = acholpre(b,A,L,R,p,Ac,B,BB,mu)

if ~exist('mu','var'), mu = 3; end

N = size(A,1);
Nc = size(Ac,1);
%% Pre-smoothing
% forward Gauss-Seidel

x = B\b;
r = b - A*x;
for k = 1:mu
    x = x + B\r;
    r = b - A*x;
end
r = r(p);     % permute into the index system of L

%% Approximated LU solve
e = L\r;
if Nc>1
    e(N-Nc+1:N) = Ac\e(N-Nc+1:N); % exact solve in the coarse grid 
end
e = R\e;      

%% Post-smoothing
e(p) = e;     % permute back to the index system of A 
x = x + e;
% backward Gauss-Seidel
for k = 1:mu+1
    r = b - A*x;
    x = x + BB\r;
end