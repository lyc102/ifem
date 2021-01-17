function pde = DataDivCurlSingularTorus2(gamma,alpha)
%%
%           curl u = g    in \Omega,
%           div (eps* u) = f    in \Omega,
%           n \cdot (eps * u) = g_N  on \Gamma_0
%  A singular vector field in a Torus like domain
%  u = curl (0,0, r^{gamma_1} sin(\alpha \theta_1)+r^{gamma_2} sin(\alpha \theta_2)) 
%  in cylindrical
%  u is in H^{gamma-\epsilon}(\Omega)

pde.Eps = @(p) [1+p(:,1).*0 0.*p(:,1) 0.*p(:,1); ...
    0*p(:,2) 1+p(:,2).*0 0*p(:,2); ...
    0*p(:,2) 0*p(:,2) 1+p(:,3).*0];
pde.exactu = @exactu;
pde.f = @(p) zeros(size(p,1),1);
pde.g = @curlu;
pde.phi = @potential;

    function s = exactu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        r2_1 = x.^2 + y.^2;
        r2_2 = (x-1).^2 + y.^2;
        theta_1 = atan2(y,x);
        theta_2 = atan2(y,x-1);
        % theta from -pi to pi, therefore the 2nd quadrant is removed in
        % Lshape as the range of alpha*theta
        
        s(:,1) = (alpha*x.*cos(alpha*theta_1)+...
            gamma(1)*y.*sin(alpha*theta_1)).*r2_1.^(-1 + gamma(1)/2)...
            + (alpha*(x-1).*cos(alpha*theta_2)+...
            gamma(2)*y.*sin(alpha*theta_2)).*r2_2.^(-1 + gamma(2)/2);
        
        s(:,2) = -(alpha*(-y).*cos(alpha*theta_1)+...
            gamma(1)*x.*sin(alpha*theta_1)).*r2_1.^(-1 + gamma(1)/2)...
            -(alpha*(-y).*cos(alpha*theta_2)+...
            gamma(2)*(x-1).*sin(alpha*theta_2)).*r2_2.^(-1 + gamma(2)/2);
        
        s(:,3) = 0*z;
    end

    function s = curlu(p) %% curl curl (0,0,phi) = -laplacian of phi
        x = p(:,1); y = p(:,2); z = p(:,3); %#ok<NASGU>
        r2_1 = x.^2 + y.^2;
        r2_2 = (x-1).^2 + y.^2;
        theta_1 = atan2(y,x);
        theta_2 = atan2(y,x-1);
        % theta from -pi to pi, therefore the 2 quadrant is removed in
        % Lshape as the range of alpha*theta
        s(:,1) = 0*x;
        s(:,2) = 0*y;
        s(:,3) = -(gamma(1)^2 - alpha^2)*r2_1.^(-1 + gamma(1)/2).*sin(alpha*theta_1)...
            -(gamma(2)^2 - alpha^2)*r2_2.^(-1 + gamma(2)/2).*sin(alpha*theta_2);
    end

    function s = potential(p)
        x = p(:,1); y = p(:,2); z = p(:,3); %#ok<NASGU>
        r2_1 = x.^2 + y.^2;
        r2_2 = (x-1).^2 + y.^2;
        theta_1 = atan2(y,x);
        theta_2 = atan2(y,x-1);
        s = r2_1.^(gamma(1)/2).*sin(alpha*theta_1) + r2_2.^(gamma(2)/2).*sin(alpha*theta_2);
    end


end