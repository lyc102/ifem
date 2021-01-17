function pde = DataDivCurlConstant
%%
%           curl u = g    in \Omega,
%           div (eps* u) = f    in \Omega,
%           n \cdot (eps * u) = g_N  on \Gamma_0
% Reference: a constant vector field

alpha = 1;

pde.Eps = @(p) [1+p(:,1)*0 p(:,3)*0 p(:,3)*0; ...
    p(:,3)*0 1+p(:,2)*0 p(:,3)*0; ...
    0*p(:,1) 0*p(:,2) 1+p(:,3)*0];


pde.exactu = @(p) alpha*[ones(size(p,1),1),...
                         ones(size(p,1),1),...
                         ones(size(p,1),1)];
pde.f = @(p) 0*p(:,1); 
pde.g = @(p) [0*p(:,1), 0*p(:,2), 0*p(:,3)];


end