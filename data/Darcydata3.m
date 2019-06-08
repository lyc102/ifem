function pde = Darcydata3
%% DARCYDATA3 Jump coefficients
%
% Use for grid with size h=1/4

pde = struct('f', @f,'g_D',@g_D,'g_N',@g_N);

% p = round(5*rand(16,1));
p = [2  3  2  3  3  4  3  5  1  2  1  0  2  1  2  4]';
pde.K = repmat(10.^(-p),2,1); % repeat for each triangle in the small square

    function s = f(pt)
    x = pt(:,1); y = pt(:,2);
    s = sin(2*pi*x).*cos(2*pi*y);
    end
    function s = g_D(pt)
    s = 0; 
    end
    function s = g_N(pt,vargin)
    s = 0; 
    end
end