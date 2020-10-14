%% DEMO of vectorization of elementwise matrix vector multiplication in MATLAB
% D(i,:,:) is a 3x3 or 2x2 matrix on the i-th element
% D can be the tensorial diffusion coefficient of a PDE
% v(i,:) is a 3x1 or 2x1 vector on the i-th element
% v can be the elementwise gradient of a certain nodal basis, e.g., 
% v(i,:) = Dphi(i,:,k) which is the k-th nodal basis's gradient
% We want to compute b(i,:) = D(i,:,:)*v(i,:)
% Even though MATLAB greatly improves the efficiency of for loop after 2016b
% using built-in vectorization operation is still faster.
% Depending on your computer the speed up is about 100x to 150x.
%
% see also Maxwell

%% set up mesh and gradient matrix
[node,elem] = cubemesh([0,1,0,1,0,1],1/32);
center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
    node(elem(:,3),:) + node(elem(:,4),:))/4;
Dlambda = gradbasis3(node,elem);
fprintf('\nSize of gradient of lambda (%d, %d, %d)\n', size(Dlambda));
NT = size(elem,1);

%% set up coefficient matrix (unfortunately arrayfun is not vectorization)
tic;
c = @(p) [1+p(:,1).^2/10 0.*p(:,1) 0.*p(:,1); ...
    0*p(:,2) 1+p(:,2).^2/10 0*p(:,2); ...
    0*p(:,2) 0*p(:,2) 1+p(:,3).^2/10];
tmp = arrayfun(@(rowidx) c(center(rowidx,:)), 1:size(center,1), 'UniformOutput',0);
c2elem = cat(3,tmp{:});
c2elem = permute(c2elem,[3,1,2]); 
% permute switch the element idx to the 1st dim
% if the Dlambda tensor has D(:,:,j) corresponds to the j-th element
% then this operation is not neccessary
fprintf('Time to generate elementwise coefficient matrix %5.4g s\n',toc);

%% benchmark the elementwise product
tic;
b = zeros(NT,3,4);
for j = 1:4
    for i = 1:size(elem,1)
        b(i,:,j) = Dlambda(i,:,j)*squeeze(c2elem(i,:,:));
    end
end
t1 = toc;
fprintf('Time to perform elementwise operation %5.4g s\n',t1);


%% benchmark the vectorized routine
tic;
b = zeros(NT,3,4);
for j = 1:4
    b(:,:,j) = sum(bsxfun(@times, c2elem, Dlambda(:,:,j)), 2);
end
t2= toc;
fprintf('Time to perform vectorized operation %5.4g s\n',t2);
fprintf('Speed up factor is %4d.\n', floor(t1/t2));