function pde = biharmonicdata
%% BIHARMONICDATA trigonometric data for biharmonic equation
%
%     f = 64*pi^4*sin(2*pi*x)*cos(2*pi*y);
%     u = sin(2*pi*x)*cos(2*pi*y);
%     Du = (2*pi*cos(2*pi*x)*cos(2*pi*y), -2*pi*sin(2*pi*x)*sin(2*pi*y));
%     w = Delta u = -8*pi^2*sin(2*pi*x)*cos(2*pi*y);
%     du/dn = g_N on [0,1]^2.
%
% Added by Jie Zhou.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'Du',@Du,'exactw',@exactw,'Dw',@Dw,...
             'g_D',@g_D,'g_N',@g_N);

% load data (right hand side function)
function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs =  64*pi^4*sin(2*pi*x).*cos(2*pi*y);
end

% exact solution
function u =  exactu(p)
    x = p(:,1); y = p(:,2);
    u =  sin(2*pi*x).*cos(2*pi*y);
end

% derivative of the exact solution
function uprime =  Du(p)
    x = p(:,1); y = p(:,2);
    uprime(:,1) = 2*pi*cos(2*pi*x).*cos(2*pi*y);
    uprime(:,2) = -2*pi*sin(2*pi*x).*sin(2*pi*y);
end

% exact solution of w = -Delta u
function w =  exactw(p)
    x = p(:,1); y = p(:,2);
    w = -8*pi^2*sin(2*pi*x).*cos(2*pi*y);
end

% derivative of the exact solution w
function wprime =  Dw(p)
    x = p(:,1); y = p(:,2);
    wprime(:,1) = -16*pi^3*cos(2*pi*x).*cos(2*pi*y);
    wprime(:,2) = 16*pi^3*sin(2*pi*x).*sin(2*pi*y);
end

% Dirichlet boundary condition
function u = g_D(p)
    u = exactu(p);
end

% Neumann boundary condition 
function f = g_N(p)
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2);
    uprime = [2*pi*cos(2*pi*p(:,1)).*cos(2*pi*p(:,2)) -2*pi*sin(2*pi*p(:,1)).*sin(2*pi*p(:,2))];
    leftbd = (abs(x)<eps);  % n = (-1,0); 
    f(leftbd) = - uprime(leftbd,1);
    rightbd = (abs(x-1)<eps); % n = (1,0); 
    f(rightbd) = uprime(rightbd,1);
    topbd = (abs(y-1)<eps);   % n = (0,1)
    f(topbd) = uprime(topbd,2);
    bottombd = (abs(y)<eps);% n = (0,-1)
    f(bottombd) = - uprime(bottombd,2);
end
end