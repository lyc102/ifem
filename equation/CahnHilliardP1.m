function [w,u,tn] = CahnHilliardP1(node,elem,pde,w0,u0,tao,times)

N = size(node,1); NT = size(elem,1);
dquadorder = 5; % the default numerical quadrature

%% Compute geometric quantities and gradient of local basis
area = simplexvolume(node,elem);
[lambda,weight] = quadpts(dquadorder);
phi = lambda;                 
nQuad = size(lambda,1);

%% Assemble stiffness and fmass matrix and jacobi matrix M0;
[A,M] = assemblematrix(node,elem);
M0 = sparse(N,N);
% assemble the Jacobi mass matrix
for i = 1:3
    for j = i:3
         M0ij=0;         
         for p = 1:nQuad 
             M0ij = M0ij + weight(p)*phi(p,i)*phi(p,j)*pde.df2du2(...
                   (phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3))));
         end 
         M0ij = M0ij.*area;                
         if (j==i) 
            M0 = M0  + sparse(elem(:,i),elem(:,j),M0ij,N,N);
         else
            M0 = M0 + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                             [M0ij; M0ij],N,N);
         end        
    end
end
% right hand side
bt = zeros(NT,3);
for p = 1:nQuad
	% quadrature points in the x-y coordinate
    for i = 1:3
        bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*pde.dfdu((...
        phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3))));
    end
end
bt = bt.*repmat(area,1,3);
rf2 = accumarray(elem(:),bt(:),[N 1]);  % the nonlinear item;
clear M0ij Dphi


global efsilonsquare
% form the linear differential operator
L = [tao*A M; M -efsilonsquare*A];   % the unknow U = [w; u]; 
rf1 = M*u0; 
rf = [rf1; rf2];       % the right hand,including two items;
U0   = [w0;u0];        % keep the two variables.
maxNewtonIt = 6;       % the max of newton iteration;
res = L*U0-rf;         % the residual ;
resnorm0 = norm(res);
% Time evolution
for k  = 1:times   
    tn = tao*k;
    %% ************************* Newton iteration *************************
    ii = 0; 
    resnorm = 1;
    while (ii < maxNewtonIt) && (resnorm > 1e-6)
       ii = ii+1;
       Jacobi = L - [sparse(N,N) sparse(N,N); sparse(N,N) M0];
       U0 = U0 - Jacobi\res;
       u0 = U0(N+1:2*N);
       M0 = sparse(N,N);  % update Jacobi matrix;
       for i = 1:3
           for j = i:3
               M0ij=0;         
               for p = 1:nQuad 
                   M0ij = M0ij + weight(p)*phi(p,i)*phi(p,j)*pde.df2du2(...
                          (phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3))));
               end
               M0ij = M0ij.*area;                
               if (j==i) 
                  M0 = M0 + sparse(elem(:,i),elem(:,j),M0ij,N,N);
               else
                  M0 = M0 + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                                   [M0ij; M0ij],N,N);
               end        
            end
       end
       % right hand side
       bt = zeros(NT,3);
       for p = 1:nQuad
           % quadrature points in the x-y coordinate
           for i = 1:3
               bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*pde.dfdu((...
                phi(p,1)*u0(elem(:,1))+phi(p,2)*u0(elem(:,2))+phi(p,3)*u0(elem(:,3))));
           end
       end
       bt = bt.*repmat(area,1,3);
       rf2 = accumarray(elem(:),bt(:),[N 1]);  % the nonlinear item;
       res = L*U0 -[rf(1:N);rf2];
       resnorm = norm(res)/resnorm0; % relative residual
    end    
    % %*********************** Newton iteration end ************************
    %% Next time step
    res(1:N) = res(1:N) + rf(1:N) - M*u0; 
    rf(1:N) = M*u0; % update right hand side for the next step
    resnorm0 = norm(res);
    %% plot solution every 5 time iterations
    disp(k)
    if mod(k,5) == 1 
        showsolution(node,elem,u0); 
        axis equal; view(2);
        pause(0.1);
    end
end
u = U0(N+1:2*N);
w = U0(1:N);
end % end of CahnHilliard
