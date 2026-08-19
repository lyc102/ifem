function pde = dataMaxwellKellogg(varargin)
%
%       curl(curl u) + beta u = J    in \Omega,
%                         n \times u = n \times g_D  on \Gamma_D,
%           n \times (mu^(-1)curl u) = n \times g_N  on \Gamma_N
%
%   Test problem for interfaces modifying the Kellogg interfaces
%   u = grad (phi)
%   where phi is the solution of the Kellogg elliptic interface problem
%   Usage: pde = dataMaxwellKellogg(gamma), gamma in {0.1, 0.5, 0.02}.
%   
%   Reference: Cai, Zhiqiang, and Shuhao Cao. "A recovery-based a posteriori error estimator for H(curl) interface problems." Computer Methods in Applied Mechanics and Engineering 296 (2015): 169-195.
%
%  Copyright (C) Long Chen. See COPYRIGHT.txt for details.



pde.mu = 1;
pde.beta = @beta;
pde.J = @J;
pde.exactu = @exactu;
pde.g_D = @exactu;
pde.curlu = [0 0 0];
pde.g_N = @g_N;
pde.divJ = 0;
pde.phi = @phi;

rho = pi/4;

if isempty(varargin) || varargin{1} == 0.1
    gamma = 0.1;
    sigma = -14.92256510455152;
    R = 161.4476387975881;
elseif varargin{1} == 0.5
    gamma = 0.5;
    sigma = -2.3561944901923448;
    R = 5.8284271247461907;
elseif varargin{1} == 0.02
    gamma = 0.02;
    sigma = -77.754418176347386;
    R = 4052.1806954768103;
else
    error('dataMaxwellKellogg:unsupportedGamma', ...
        'Unsupported gamma = %g. Supported values are 0.1, 0.5, and 0.02.', varargin{1});
end

    function s = beta(p)
        idx = (p(:,1).*p(:,2) >0);
        s  = ones(size(p,1),1);
        s(idx) = R;
    end

    function s = g_N(p, n)
        s = zeros(size(p,1),3);
    end

    function s = exactu(p)
        % the gradient of exact solution phi for Kellogg
        theta = atan2(p(:,2),p(:,1));
        theta = (theta>= 0).*theta +(theta<0).*(theta+2*pi);
        t = 1+(p(:,2).^2)./(p(:,1).^2);
        r = sqrt(p(:,2).^2+p(:,1).^2);
        rg = r.^gamma;
        
        ux1 = (p(:,1)> 0.0 & p(:,2)>= 0.0).*...
            (rg.*gamma./r.*cos((pi/2-sigma)*gamma)./r.*p(:,1)...
            .*cos((theta-pi/2+rho)*gamma)...
            +rg.*cos((pi/2-sigma)*gamma).*sin((theta-pi/2+rho)*gamma)...
            *gamma.*p(:,2)./(p(:,1).^2)./t);
        
        uy1 = (p(:,1)> 0.0 & p(:,2)>= 0.0).*...
            (rg*gamma./r.*cos((pi/2-sigma).*gamma)...
            .*cos((theta-pi/2+rho).*gamma)./r.*p(:,2)...
            -rg.*cos((pi/2-sigma).*gamma).*sin((theta-pi/2+rho)*gamma)...
            *gamma./p(:,1)./t);
        
        ux2 = (p(:,1)<= 0.0 & p(:,2)> 0.0).*...
            (rg*gamma./r.*cos(rho.*gamma)./r.*p(:,1)...
            .*cos((theta-pi+sigma).*gamma)...
            +rg.*cos(rho*gamma).*sin((theta-pi+sigma).*gamma)*gamma...
            .*p(:,2)./(p(:,1).^2)./t);
        
        uy2 = (p(:,1)<= 0.0 & p(:,2)> 0.0).*...
            (rg*gamma./r.*cos(rho.*gamma)./r.*p(:,2)...
            .*cos((theta-pi+sigma)*gamma)-rg.*cos(rho*gamma)...
            .*sin((theta-pi+sigma).*gamma)*gamma./p(:,1)./t);
        
        ux3 = (p(:,1)< 0.0 & p(:,2)<= 0.0).*...
            (rg*gamma./r.*cos(sigma.*gamma)./r.*p(:,1)...
            .*cos((theta-pi-rho).*gamma) ...
            +rg.*cos(sigma.*gamma).*sin((theta-pi-rho).*gamma).*gamma...
            .*p(:,2)./(p(:,1).^2)./t);
        
        uy3 = (p(:,1)< 0.0 & p(:,2)<= 0.0).*...
            (rg*gamma./r.*cos(sigma.*gamma)./r.*p(:,2)...
            .*cos((theta-pi-rho).*gamma) -rg.*cos(sigma*gamma)...
            .*sin((theta-pi-rho).*gamma)*gamma./p(:,1)./t);
        
        ux4 = (p(:,1)>= 0.0& p(:,2)< 0.0).*...
            (rg*gamma./r.*cos((pi/2-rho).*gamma)./r.*p(:,1)...
            .*cos((theta-3*pi/2-sigma)*gamma) ...
            +rg.*cos((pi/2-rho).*gamma).*sin((theta-3*pi/2-sigma)*gamma)...
            *gamma.*p(:,2)./(p(:,1).^2)./t);
        
        uy4 = (p(:,1)>= 0.0 & p(:,2)< 0.0).*...
            (rg*gamma./r.*cos((pi/2-rho).*gamma)./r.*p(:,2)...
            .*cos((theta-3*pi/2-sigma)*gamma)...
            -rg.*cos((pi/2-rho)*gamma).*sin((theta-3*pi/2-sigma)*gamma)...
            *gamma./p(:,1)./t);
        
        s(:,1) = ux1+ux2+ux3+ux4;
        s(:,2) = uy1+uy2+uy3+uy4;
        s(:,3) = 0*p(:,3);
    end

    function s = J(p)
        s = repmat(beta(p),1,3).*exactu(p);
    end

    function s = phi(p)
        r = sqrt(p(:,2).^2+p(:,1).^2);
        theta = atan2(p(:,2),p(:,1));
        theta = (theta>= 0).*theta + (theta<0).*(theta+2*pi);
        mu = (theta>= 0 & theta<pi/2).*cos((pi/2-sigma)*gamma).*cos((theta-pi/2+rho)*gamma)...
            +(theta>= pi/2 & theta<pi).*cos(rho*gamma).*cos((theta-pi+sigma)*gamma)...
            +(theta>= pi & theta<1.5*pi).*cos(sigma*gamma).*cos((theta-pi-rho)*gamma)...
            +(theta>= 1.5*pi & theta<2*pi).*cos((pi/2-rho)*gamma).*cos((theta-1.5*pi-sigma)*gamma);
        s =  r.^gamma.*mu;
    end

end
