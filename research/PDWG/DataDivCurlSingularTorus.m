function pde = DataDivCurlSingularTorus(gamma,alpha,varargin)
%%
%           curl u = g    in \Omega,
%           div (eps* u) = f    in \Omega,
%           n \cdot (eps * u) = g_N  on \Gamma_0
%  A singular vector field in a Torus like domain
%  u = curl (0,0, r^{gamma} sin(\alpha \theta)) in cylindrical
%  u is in H^{gamma-\epsilon}(\Omega)


if length(gamma) > 1; gamma = gamma(1); end
beta = 0;
if nargin > 2; beta = varargin{1}; end

pde.Eps = @(p) [1+p(:,1).*0 0.*p(:,1) 0.*p(:,1); ...
    0*p(:,2) 1+p(:,2).*0 0*p(:,2); ...
    0*p(:,2) 0*p(:,2) 1+p(:,3).*0];
pde.exactu = @exactu;
pde.f = @(p) beta*pi*cos(pi*p(:,1)).*cos(pi*p(:,2));
pde.g = @curlu;
pde.phi = @potential;

    function s = exactu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        r2 = x.^2 + y.^2;
        theta = atan2(y,x);
        % theta from -pi to pi, therefore the 2nd quadrant is removed in
        % Lshape as the range of alpha*theta
        
        s(:,1) = (alpha*x.*cos(alpha*theta)+...
            gamma*y.*sin(alpha*theta)).*r2.^(-1 + gamma/2);
        
        s(:,2) = -(alpha*(-y).*cos(alpha*theta)+...
            gamma*x.*sin(alpha*theta)).*r2.^(-1 + gamma/2);
        
        s(:,3) = 0*z;
        
        s(:,1) = s(:,1) + beta*cos(pi*y).*sin(pi*x);
        s(:,2) = s(:,2) - beta*cos(pi*x).*sin(pi*y);
    end

    function s = curlu(p) %% curl curl (0,0,phi) = -laplacian of phi
        x = p(:,1); y = p(:,2); z = p(:,3); %#ok<NASGU>
        r2 = x.^2 + y.^2;
        theta = atan2(y,x);
        % theta from -pi to pi, therefore the 2 quadrant is removed in
        % Lshape as the range of alpha*theta
        s(:,1) = 0*x;
        s(:,2) = 0*y;
        s(:,3) = -(gamma^2 - alpha^2)*r2.^(-1 + gamma/2).*sin(alpha*theta);
        s(:,3) = s(:,3) + 2*beta*pi*sin(pi*x).*sin(pi*y);
    end

    function s = potential(p)
        x = p(:,1); y = p(:,2); z = p(:,3); %#ok<NASGU>
        r2 = x.^2 + y.^2;
        theta = atan2(y,x);
        s = r2.^(gamma/2).*sin(alpha*theta);
    end


end