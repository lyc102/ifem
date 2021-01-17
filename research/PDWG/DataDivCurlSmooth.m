function pde = DataDivCurlSmooth(varargin)
%%
%           curl u = g    in \Omega,
%           div (eps* u) = f    in \Omega,
%           n \cdot (eps * u) = g_N  on \Gamma_0
% Reference: a simple solenoidal vector field

if isempty(varargin)
    alpha = 1/4; % strength of solenoidal part cos*sin
    beta = 1/4; % strength of irrotational part (x,y,z)
else
    alpha = varargin{1};
    beta = varargin{2};
end

pde.Eps = @(p) [3+p(:,1)*0 p(:,3)*0 p(:,3)*0; ...
    p(:,3)*0 2+p(:,2)*0 p(:,3)*0; ...
    0*p(:,1) 0*p(:,2) 1+p(:,3)*0];


pde.exactu = @(p) alpha*[cos(pi*p(:,2)).*sin(pi*p(:,1)),...
                  -cos(pi*p(:,1)).*sin(pi*p(:,2)), 0*p(:,3)] + ...
                  beta*[p(:,1), p(:,2), p(:,3)];
pde.f = @(p) alpha*pi*cos(pi*p(:,1)).*cos(pi*p(:,2)) + beta*6; 
pde.g = @(p) alpha*[0*p(:,1), 0*p(:,2), 2*pi*sin(pi*p(:,1)).*sin(pi*p(:,2))];


end