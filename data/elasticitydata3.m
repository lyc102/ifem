function pde = elasticitydata3(para)
%% ELASTICITYDATA3 data for the elasticity problem in three dimensions
%
% Modified from elasticitydata by Huayi Wei.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Lame constants
if nargin == 0
    lambda = 1;
    mu = 1;
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'lambda') || isempty(para.lambda)
        lambda = 1;
    else
        lambda = para.lambda;
    end
    if ~isfield(para,'mu') || isempty(para.mu)
        mu = 1;
    else
        mu = para.mu;
    end
end

pde = struct('lambda',lambda,'mu',mu, 'f', @f, 'exactu',@exactu,'g_D',@g_D);
%%%%%% subfunctions %%%%%%
    function s = f(p)
% f = - mu \Delta u - (mu + lambda) grad(div u) 
       x = p(:,1); y = p(:,2); z = p(:,2);
       sx = sin(pi*x); sy = sin(pi*y); sz = sin(pi*z);
       cx = cos(pi*x); cy = cos(pi*y); cz = cos(pi*z);
       f1 = (mu*4*pi^2)*exactu(p);
       f2 = zeros(size(p));
       f2(:,1) = (mu+lambda)*pi^2*(-sx.*sy.*sz.*cz - sx.*cy.*sz.*sx + ...
                                    cx.*cy.*sz.*cx + cx.*sy.*cz.*cy);
       f2(:,2) = (mu+lambda)*pi^2*(cx.*cy.*sz.*cz - sx.*sy.*sz.*cx + ...
                                   sx.*cy.*cz.*cy - sx.*sy.*cz.*sy);
       f2(:,3) = (mu+lambda)*pi^2*(cx.*sy.*cz.*cz - cx.*sy.*sz.*sz + ...
                                   sx.*cy.*cz.*cx - sx.*sy.*sz.*cy); 
       s = f1 - f2;
    end
    function s = exactu(p)
       x = p(:,1); y = p(:,2); z = p(:,2);
       temp = sin(pi*x).*sin(pi*y).*sin(pi*z); % zero boundary condition
       s = [cos(pi*z).*temp cos(pi*x).*temp cos(pi*y).*temp];    
    end
    function s = g_D(p)
       s = exactu(p); 
    end
end