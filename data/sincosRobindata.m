function pde = sincosRobindata
%% SINCOSDATA trigonometric  data for Poisson equation
%
%     f = 8*pi^2*sin(2*pi*x)*cos(2*pi*y);
%     u = sin(2*pi*x)*cos(2*pi*y);
%     Du = (2*pi*cos(2*pi*x)*cos(2*pi*y), -2*pi*sin(2*pi*x)*sin(2*pi*y));
%
%  Impose Robin boundary condition g_R*u+ du/dn = g_N on [0,1]^2.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'Du',@Du,'g_D',@g_D,'g_R',@g_R,'g_N',@g_N);

% load data (right hand side function)
function rhs =  f(p)
    x = p(:,1); y = p(:,2);
    rhs =  8*pi^2*sin(2*pi*x).*cos(2*pi*y);
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
% Dirichlet boundary condition
function u =  g_D(p)
    u =  exactu(p);
end
% coefficients of Robin boundary condition g_R*u+ du/dn = g_N 
function f = g_R(p)
    f = ones(size(p,1),1);
end
% 
function f = g_N(p)
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2);
    uprime = Du(p);
    leftbd = (abs(x)<eps);  % n = (-1,0); 
    f(leftbd) = g_R(p(leftbd,:)).*exactu(p(leftbd,:)) - uprime(leftbd,1);
    rightbd = (abs(x-1)<eps); % n = (1,0); 
    f(rightbd) = g_R(p(rightbd,:)).*exactu(p(rightbd,:)) + uprime(rightbd,1);
    topbd = (abs(y-1)<eps);   % n = (0,1)
    f(topbd) = g_R(p(topbd,:)).*exactu(p(topbd,:)) + uprime(topbd,2);
    bottombd = (abs(y)<eps);% n = (0,-1)
    f(bottombd) = g_R(p(bottombd,:)).*exactu(p(bottombd,:)) - uprime(bottombd,2);
end
end