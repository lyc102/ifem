function pde = Maxwelldata5
%% MAXWELLDATA5 linear polynomial data
%
%       curl(mu^(-1)curl u) - epsilon u = J    in \Omega,  
%                                 n × u = n × g_D  on \Gamma_D,
%                   n × (mu^(-1)curl u) = n × g_N  on \Gamma_N?
%
% mu = 1; 
% epsilon = - 1; % positive type
% - u satisfies 0 Dirichlet boundary condition at cube [-1,1]^3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

pde.mu = 1;
pde.epsilon = -1;
pde.J = @J;
pde.exactu = @exactu;
pde.g_D = @exactu;
pde.curlu = @curlu;

    function s = J(p)
    % x = p(:,1); y = p(:,2); z = p(:,3);
    % s = [-sin(y), 0*x, 0*z]; 
    s = repmat([2, 0, 0],size(p,1),1) + exactu(p);
    end

    function s = exactu(p)
    x = p(:,1); y = p(:,2); z = p(:,3);
    s = [-y.^2, x, 0*z];
    end

    function s = curlu(p)
    x = p(:,1); y = p(:,2); z = p(:,3); %#ok<NASGU>
    s = [0*x, 0*y, 1+2*y];
    % s = repmat([0, 0, 2],size(p,1),1);
    end
end