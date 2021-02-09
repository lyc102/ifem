function pde = quadCurlDataLshape2(varargin)
%%
%           curl^4 u = f    in \Omega,
%           div u = 0    in \Omega,
%           u x n and (curl u) x n are both nonhomogeneous
%  Reference: adapted from Example 3 in
%  H. Guo, Z. Zhang, Q. Zou: A C0 Linear Finite Element Method for Biharmonic Problems
%  The example 3 in this paper is wrong (either typo or error)
%  u = curl^2 (0,0,phi) = (0,0,-\Delta phi), curl^2 u = curl (0,0, \Delta^2 phi)
%  curl^4 u = (0,0, -Delta^3 phi) 
%  phi = r^{\gamma} \sin(\alpha * \theta)
%  (-1,1)^2\times (0,1/2)\ (x>0 and y<0)
%  curl^2 u = 0 if \gamma = \alpha+2 when r >0
%  
%  Reference: example 3 in Error analysis of a decoupled finite element method for quad-curl problems
%  https://arxiv.org/abs/2102.03396
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.


pde = struct('exactu', @exactu,'g_D', @g_D, 'g', @divu, ...
    'curlu', @curlu, 'curlu_D',@curlu_D, 'curlcurlu',@curlcurlu, ...
    'tricurlu',@tricurlu,'quadcurlu',@quadcurlu);

if isempty(varargin)
    gamma = 8/3;
    alpha = 2/3;
else
    gamma = varargin{1};
    alpha = varargin{2};
end

    function s = exactu(p)
        x = p(:,1); y = p(:,2);
        s(:,1) = zeros(size(p,1),1);
        s(:,2) = s(:,1);
        r2 = x.^2 + y.^2;
        theta = atan2(y,x);
        theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        s(:,3) = (alpha^2 - gamma^2)*r2.^(gamma/2-1).*sin(alpha*theta);
    end

    function s = curlu(p)
        x = p(:,1); y = p(:,2);
        r2 = x.^2 + y.^2;
        theta = atan2(y,x);
        theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        s(:,1) = (alpha^2 - gamma^2)*(alpha*x.*cos(alpha*theta)+...
            (gamma-2)*y.*sin(alpha*theta)).*r2.^(gamma/2-2);
        
        s(:,2) = (alpha^2 - gamma^2)*(-(gamma-2)*x.*sin(alpha*theta)+...
            alpha*y.*cos(alpha*theta)).*r2.^(gamma/2-2);
        
        s(:,3) = zeros(size(p,1),1);
    end

    function s = curlcurlu(p)
        x = p(:,1); y = p(:,2); 
        s(:,1) = zeros(size(p,1),1);
        s(:,2) = s(:,1);
        r2 = x.^2 + y.^2 + 2*eps;
        theta = atan2(y,x);
        theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        s(:,3) = (alpha^2 - gamma^2)*(alpha^2 - (gamma-2)^2)...
            *r2.^(gamma/2-2).*sin(alpha*theta);
    end


    function s = tricurlu(p)
        x = p(:,1); y = p(:,2);
        r2 = x.^2 + y.^2 + 2*eps;
        theta = atan2(y,x);
        theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        
        s(:,1) = (alpha^2 - gamma^2)*(alpha^2 - (gamma-2)^2)*...
            (alpha*x.*cos(alpha*theta)+(gamma-4)*y.*sin(alpha*theta)).*r2.^(gamma/2-3);
        
        s(:,2) = (alpha^2 - gamma^2)*(alpha^2 - (gamma-2)^2)*...
            (-(gamma-4)*x.*sin(alpha*theta)+alpha*y.*cos(alpha*theta)).*r2.^(gamma/2-3);
       
        s(:,3) = zeros(size(p,1),1);
    end

    function s = quadcurlu(p)
        x = p(:,1); y = p(:,2);
        s(:,1) = zeros(size(p,1),1);
        s(:,2) = s(:,1);
        r2 = x.^2 + y.^2 + 2*eps;
        theta = atan2(y,x);
        theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
        s(:,3) = (alpha^2 - gamma^2)*(alpha^2 - (gamma-2)^2)*(alpha^2 - (gamma-4)^2)...
            *r2.^(gamma/2-3).*sin(alpha*theta);
    end


    function s = g_D(p) % Dirichlet boundary condition for u
        s = exactu(p);
    end

    function s = divu(p)
        s = zeros(size(p,1),1);
    end
end
