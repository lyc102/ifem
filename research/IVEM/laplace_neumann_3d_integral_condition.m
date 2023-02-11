% clear;
% clc;


% Solve this linear equations:
%   Au = fh^2 + 2gh


%% global parameters
x1 = -1; x2 = 1;
y1 = -1; y2 = 1;
z1 = -1; z2 = 1;
M = 81;
N = M;
L = M;
h = (x2-x1) / (N-1);
k = (y2-y1) / (M-1);
l = (z2-z1) / (L-1);
x = x1:h:x2;
y = y1:k:y2;
z = z1:l:z2;
% z = zeros(L);

mu = 0;


%% Exact solution
ue = zeros(M, N, L);
for i = 1:M
    ue(i, :, :) = sin(x(i) + y + z');
end


%% Construct the coefficient matrix
% Build matrix B
r = 1/3*(6 + mu*h^2) * ones(M, 1);
% r = 0.5*4 * ones(M, 1);
r1 = -1*ones(M-1, 1);
B = diag(r, 0) + diag(r1, 1) + diag(r1, -1);
B(1, 2) = -2;
B(M, M-1) = B(1, 2);

% Sparse matrix B
B = sparse(B);

% Build sparse identity matrix
I = eye(M);
I = sparse(I);

% I2 = eye(M*M);
% I2 = sparse(I2);

% Build tridiagonal(five diagonal) block matrix A
% A = kron(B, I) + kron(I, B);

% Build seven-diagonal matrix
% C = kron(A, I) + kron(I2, B);
C = kron(I, kron(I, B)) + kron(I, kron(B, I)) + kron(kron(B, I), I);


%% Get f and g
% Get f
% FF1 = @(x, y) (mu(1)+2) * sin(x + y);
% FF2 = @(x, y) (mu(2)+2) * sin(x + y);

FF = @(x, y, z) (mu+3) * sin(x+y+z);
F = zeros(M, N, L);
for i = 1:M
    F(i, :, :) = FF(x(i), y, z');
end

% update A and F
% for i = 1:N
%     for j = 1:M
%         if board(i, j) == 1
% %             disp("****")
% %             (j-1)*(N+1) + i
% %             (x(j)-c(1,1))^2/(a(1)^2) + (y(i)-c(1,2))^2/(b(1)^2)
% %             r(1)^2
%             A((j-1)*N + i, (j-1)*N + i) = 4 + mu(2)*h^2;
% %             F(i, j) = FF2(x(j), y(i));
%         else
% %             F(i, j) = FF1(x(j), y(i));
%             continue;
%         end
%     end
% end

% F
f = reshape(F, M*N*L, 1);

% Get g
G = zeros(M, N, L);

G(1, :, :) = -cos(y+z'-1);
G(M, :, :) = cos(y+z'+1);
G(:, 1, :) = -cos(x+z'-1);
G(:, N, :) = cos(x+z'+1);
G(:, :, 1) = -cos(x+y'-1);
G(:, :, L) = cos(x+y'+1);

G(1, 1, :) = -cos(z-2) - cos(z-2);
G(1, N, :) = -cos(z) + cos(z);
G(M, 1, :) = cos(z) - cos(z);
G(M, N, :) = cos(z+2) + cos(z+2);
G(1, :, 1) = -cos(y-2) - cos(y-2);
G(1, :, L) = -cos(y) + cos(y);
G(M, :, 1) = cos(y) - cos(y);
G(M, :, L) = cos(y+2) + cos(y+2);
G(:, 1, 1) = -cos(x-2) - cos(x-2);
G(:, 1, L) = -cos(x) + cos(x);
G(:, N, 1) = cos(x) - cos(x);
G(:, N, L) = cos(x+2) + cos(x+2);

G(1, 1, 1) = -cos(-3) - cos(-3) - cos(-3);
G(1, N, 1) = -cos(-1) + cos(-1) - cos(-1);
G(M, 1, 1) = cos(-1) - cos(-1) - cos(-1);
G(1, 1, L) = -cos(-1) - cos(-1) + cos(-1);
G(M, N, 1) = cos(1) + cos(1) - cos(1);
G(M, 1, L) = cos(1) - cos(1) + cos(1);
G(1, N, L) = -cos(1) + cos(1) + cos(1);
G(M, N, L) = cos(3) + cos(3) + cos(3);

% G(1, 1, :) = G(1, 2, :) + G(2, 1, :);
% G(1, N, :) = G(1, N-1, 1) + G(2, N, 1);
% G(M, 1, :) = G(M-1, 1, 1) + G(M, 2, 1);
% G(M, N, :) = G(M-1, N, 1) + G(M, N-1, 1);
% G(1, :, 1) = G(2, :, 1) + G(1, :, 2);
% G(1, :, L) = G(2, :, L) + G(1, :, L-1);
% G(M, :, 1) = G(M-1, :, 1) + G(M, :, 2);
% G(M, :, L) = G(M-1, :, L) + G(M, :, L-1);
% G(:, 1, 1) = G(:, 2, 1) + G(:, 1, 2);
% G(:, 1, L) = G(:, 2, L) + G(:, 1, L-1);
% G(:, N, 1) = G(:, N-1, 1) + G(:, N, 2);
% G(:, N, L) = G(:, N-1, L) + G(:, N, L-1);


% G(1, 1, 1) = 3*G(1, 1, 1);
% G(1, N, 1) = 3*G(1, N, 1);
% G(M, 1, 1) = 3*G(M, 1, 1);
% G(1, 1, L) = 3*G(1, 1, L);
% G(M, N, 1) = 3*G(M, N, 1);
% G(M, 1, L) = 3*G(M, 1, L);
% G(1, N, L) = 3*G(1, N, L);
% G(M, N, L) = 3*G(M, N, L);



g = reshape(G, M*N*L, 1);


%% integration condition
tau = zeros(M, N, L);
idx1 = 1:2:N;
idx2 = 2:2:N;
% tau(1, idx1) = h/3;
% tau(1, idx2) = 4*h/3;
% tau(M, :) = tau(1, :);
% tau(:, 1) = tau(1, :);
% tau(:, N) = tau(1, :);

tau(idx1, idx1, 1) = h^2/9;
tau(idx1, idx2, 1) = 4*h^2/9;
tau(idx2, idx1, 1) = 4*h^2/9;
tau(idx2, idx2, 1) = 16*h^2/9;
tau(:, :, L) = tau(:, :, 1);
tau(:, 1, :) = tau(:, :, 1);
tau(:, N, :) = tau(:, :, 1);
tau(1, :, :) = tau(:, :, 1);
tau(M, :, :) = tau(:, :, 1);

ssize = size(f, 1);
tau_vec = reshape(tau, M*N*L, 1);
% A(ssize+1, :) = tau_vec;
C(ssize+1, :) = tau_vec';
C(:, ssize+1) = [tau_vec; 1];

rhs = (h^2 * f + 2.0 * h * g);
rhs(ssize+1, 1) = 0;


%% Compute u
fprintf("Solve\n");
tic
% u = A \ (h^2 * f + 2 * h * g);

% u = C \ rhs;

% LL = ichol(C);
alpha = 5;
LL = ichol(C,struct('type','ict','droptol',1e-3,'diagcomp',alpha));
[u, flag, relres, iter, resvec] = pcg(C, rhs, 1e-8, 1000, LL, LL');


toc
u = u(1:end-1);
U = reshape(u, M, N, L);

% f = zeros(N, M);
% f(1, :) = U(1, :);
% f(N, :) = U(N, :);
% f(:, 1) = U(:, 1);
% f(:, M) = U(:, M);
% 
% f = sparse(f);


%% Compute error
error_inf_norm = max(max(max(abs(U - ue))));
error_l2_norm = sqrt(sum(sum(sum((U - ue).^2))) / M / N / L);
% error = sqrt(sum(sum((U - ue).^2)));

fprintf("Inf-norm: %f\n", error_inf_norm);
fprintf("L2-norm: %f\n", error_l2_norm);

% order
% order1 = log(error_20/error_40) / log(h_20/h_40);


% Plotting the solution
% figure()
% surf(x, y, U, 'FaceAlpha', 0.5);
% surf(x, y, U, 'EdgeColor', 'none');
% hold on
% surf(x, y, ue,'EdgeColor','none');
% % surf(x,y,p - u,'EdgeColor','none');
% shading interp
% title({'2-D Laplace''s equation'})
% xlabel('Spatial co-ordinate (x) \rightarrow')
% ylabel('{\leftarrow} Spatial co-ordinate (y)')
% zlabel('Solution profile (P) \rightarrow')


%% return
% Un = U;
% coef_mat = A;

