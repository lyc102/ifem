function pde = sincosdata3Robin1
%% SINCOSDATA3ROBIN1 trigonometric data for Poisson equation in 3-D
%
%     f = 3*sin(x)*sin(y)*sin(z);
%     u = sin(x)*sin(y)*sin(z);
%     Du = (cos(x)*sin(y)*sin(z), 
%           sin(x)*cos(y)*sin(z),
%           sin(x)*sin(y)*cos(z));
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde = struct('f',@f,'exactu',@exactu,'g_D',@exactu,'g_R',@g_R,'Du',@Du,'g_N',@g_N);
eps = 1e-10;
    % load data (right hand side function)
function s = f(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = 3*sin(x).*sin(y).*sin(z);
end
    % exact solution
function s = exactu(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = sin(x).*sin(y).*sin(z); 
end

    % Derivative of the exact solution
function s = Du(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s(:,1) = cos(x).*sin(y).*sin(z);
    s(:,2) = sin(x).*cos(y).*sin(z);
    s(:,3) = sin(x).*sin(y).*cos(z);
end

% coefficients of Robin boundary condition g_R*u+ du/dn = g_N 
function s = g_R(p)
    s = ones(size(p,1),1);
end
%
function f = g_N(p)
    f = zeros(size(p,1),1);
    x = p(:,1); y = p(:,2); z = p(:,3);
    uprime = [cos(x).*sin(y).*sin(z) ...
              sin(x).*cos(y).*sin(z) ...
              sin(x).*sin(y).*cos(z)];
    downbd = (abs(z)<eps); % n = (0,0,-1)
    f(downbd) = g_R(p(downbd,:)).*exactu(p(downbd,:)) - uprime(downbd,3);
    upbd = (abs(z-pi)<eps);% n = (0,0,1)
    f(upbd) = g_R(p(upbd,:)).*exactu(p(upbd,:)) + uprime(upbd,3);
    leftbd = (abs(x)<eps); % n = (-1,0,0)
    f(leftbd) = g_R(p(leftbd,:)).*exactu(p(leftbd,:)) - uprime(leftbd,1);
    rightbd = (abs(x-pi)<eps); % n = (1,0,0)
    f(rightbd) = g_R(p(rightbd,:)).*exactu(p(rightbd,:)) + uprime(rightbd,1);
    frontbd = (abs(y)<eps); % n = (0,-1,0)
    f(frontbd) = g_R(p(frontbd,:)).*exactu(p(frontbd,:)) - uprime(frontbd,2);
    backbd = (abs(y-pi)<eps); % n =(0,1,0)
    f(backbd) = g_R(p(backbd,:)).*exactu(p(backbd,:)) + uprime(backbd,2);
end
end