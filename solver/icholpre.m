function x = icholpre(b,A,L,R)

%% Pre-smoothing by Gauss-Seidel
x = tril(A)\b;
r = b - A*x;

%%
e = L\r;
e = R\e;      % solve by achol decomposition L*L'
x = x + e;

%% Post-smoothing
x = x + triu(A)\(b-A*x);