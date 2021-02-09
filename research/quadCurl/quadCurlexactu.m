%% A script to generate copy'able functions for quad-curl problem
% used in "Error analysis of a decoupled finite element method for quad-curl problems"
% https://arxiv.org/abs/2102.03396
% 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear; clc;
syms x y z
% U = [sin(y), -sin(z), sin(x)];
% U = [0, 0, (sin(x))^2*(sin(y))^2*sin(z)];
U = [0, 0, (sin(x)).^3.*(sin(y)).^3.*(sin(z)).^2];
% U = [0, 0, x^2*(1-x)^2*y^2*(1-y)^2*z^2*(1-z)^2];

X = [x y z];

%%
curlu = curl(U,X);
curlcurlu = curl(curlu,X);
tricurlu = curl(curlcurlu,X);
quadcurlu = curl(tricurlu,X);
pentacurlu = curl(quadcurlu,X);
hexacurlu = curl(pentacurlu,X);

vectorsStr = {'curlu', 'curlcurlu', 'tricurlu', ...
              'quadcurlu', 'pentacurlu', 'hexacurlu'};

%%
for idx = 1:length(vectorsStr)
    vector = eval(vectorsStr{idx});
    fprintf("function s = %s(p)\n",vectorsStr{idx})
    fprintf("x = p(:,1); y = p(:,2); z = p(:,3);\n")
    for j = 1:3
        vStr = char(vector(j));
        vStr = strrep(vStr, '*', '.*');
        vStr = strrep(vStr, '^', '.^');
        fprintf("s(:,%d) = %s;\n", j, vStr)
    end
    fprintf('end\n\n')
end