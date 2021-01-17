function pde = DataDivCurlLinear(varargin)
%%
%           curl u = g    in \Omega,
%           div (eps* u) = f    in \Omega,
%           n \cdot (eps * u) = g_N  on \Gamma_0
% Reference: a linear vector field with modifiable curl and div

if isempty(varargin)
    alpha = 1; % strength of irrotational part (x,y,z)
    beta = 1; % strength of solenoidal part (z,x,y)
else
    alpha = varargin{1};
    beta = varargin{2};
end

pde.Eps = @(p) [1+p(:,1)*0 p(:,3)*0 p(:,3)*0; ...
    p(:,3)*0 1+p(:,2)*0 p(:,3)*0; ...
    0*p(:,1) 0*p(:,2) 1+p(:,3)*0];


pde.exactu = @(p) alpha*[p(:,1), p(:,2), p(:,3)] + ...
                   beta*[p(:,3), p(:,1), p(:,2)];
pde.f = @(p) 3*alpha*ones(size(p,1),1); 
pde.g = @(p) beta*[ones(size(p,1),1),...
                  ones(size(p,1),1),...
                  ones(size(p,1),1)];


end