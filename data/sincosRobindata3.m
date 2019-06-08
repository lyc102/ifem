function pde = sincosRobindata3
%% SINCOSDATA3 trigonometric data for Poisson equation in 3-D
%
%     f = 3*pi^2*cos(pi*x)*cos(pi*y)*cos(pi*z);
%     u = cos(pi*x)*cos(pi*y)*cos(pi*z);
%     Du = (-pi*sin(pi*x)*cos(pi*y)*cos(pi*z), 
%           -pi*cos(pi*x)*sin(pi*y)*cos(pi*z),
%           -pi*cos(pi*x)*cos(pi*y)*sin(pi*z));
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'g_R',@g_R,'Du',@Du,'g_N',@g_N);

    % load data (right hand side function)
function s = f(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = 3*pi^2*cos(pi*x).*cos(pi*y).*cos(pi*z);
end
    % exact solution
function s = exactu(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = cos(pi*x).*cos(pi*y).*cos(pi*z); % for neumann boundary condition, int_u =0
end

    % Derivative of the exact solution
function s = Du(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s(:,1) = -pi*sin(pi*x).*cos(pi*y).*cos(pi*z);
    s(:,2) = -pi*cos(pi*x).*sin(pi*y).*cos(pi*z);
    s(:,3) = -pi*cos(pi*x).*cos(pi*y).*sin(pi*z);
end

    % Dirichlet boundary condition
function s = g_D(p)
    s = exactu(p);
end

    % coefficients of Robin boundary condition g_R*u+ du/dn = g_N 
function s = g_R(p)
    s = ones(size(p,1),1);
end
%
function f = g_N(p)
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2); z = p(:,3);
    uprime = [-pi*sin(pi*x).*cos(pi*y).*cos(pi*z) ...
              -pi*cos(pi*x).*sin(pi*y).*cos(pi*z) ...
              -pi*cos(pi*x).*cos(pi*y).*sin(pi*z)];
    downbd = (abs(z)<eps); % n = (0,0,-1)
    f(downbd) = g_R(p(downbd,:)).*exactu(p(downbd,:)) - uprime(downbd,3);
    upbd = (abs(z-1)<eps);% n = (0,0,1)
    f(upbd) = g_R(p(upbd,:)).*exactu(p(upbd,:)) + uprime(upbd,3);
    leftbd = (abs(x)<eps); % n = (-1,0,0)
    f(leftbd) = g_R(p(leftbd,:)).*exactu(p(leftbd,:)) - uprime(leftbd,1);
    rightbd = (abs(x-1)<eps); % n = (1,0,0)
    f(rightbd) = g_R(p(rightbd,:)).*exactu(p(rightbd,:)) + uprime(rightbd,1);
    backbd = (abs(y-1)<eps); % n =(0,1,0)
    f(backbd) = g_R(p(backbd,:)).*exactu(p(backbd,:)) + uprime(backbd,2);
    frontbd = (abs(y)<eps); % n = (0,-1,0)
    f(frontbd) = g_R(p(frontbd,:)).*exactu(p(frontbd,:)) - uprime(frontbd,2);
end
end