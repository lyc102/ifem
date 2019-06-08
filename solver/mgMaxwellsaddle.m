function [u,p,info,info2] =  mgMaxwellsaddle(A,G,f,g,node,elem,bdFlag,Me,grad,option)
%Solve the maxwell system with divgence free condition,
%         [A  G] [u]  = f                (1.1)
%         [G' O] [p]  = g               (1.2)
% where  G = M_e*grad.
% This system can be rewritten as 
%         [A+G*DMinv*G'       G] [u]  = f +G*DMinv*g
%         [G'                 O] [p]  = g
% in fact that
%  [A+G*DMinv*G'  G] [I    grad]                = [Abar   O  ]
%  [G'            O] [0  -Dminv*grad'*M_e*grad ]= [G'     A_p]
%
%  We solve the (1.1)-(1.2) with the GMRES, the preconditioner is
%
%       [I  grad                ]  [Abar O   ]^{-1}
%       [O  -Dminv*grad'*M_e*grad] [G'   A_p]
%
% Created by Long chen and Jie Zhou on Aug,2015.


Nf = length(f); 
Ng = length(g);
t = cputime;

%% Parameters
if ~exist('option','var'), option = []; end
Ndof = Nf + Ng; 
option = mgoptions(option,Ndof);    % parameters
d = size(node,2);

%% Set up auxiliary matrices
if d == 2
    area = simplexvolume(node,elem);
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,...
                        [max(elem(:)),1]);
elseif d == 3
    volume = abs(simplexvolume(node,elem)); % uniform refine in 3D is not orientation presereved
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volume;volume;volume;volume]/4,[max(elem(:)),1]);    
end
DMinv = spdiags(1./Mvlump(option.isFreeNode),0,Ng,Ng);
f = f + G*(DMinv*g);  % add second equation to the first one
Abar = A + G*DMinv*G'; % Hodge Laplacian
Ap = grad'*Me*grad;    % scalar Laplacian

%% Solve 
option.printlevel = 1;
u0 = option.x0(1:Nf);
p0 = option.x0(Nf+1:end);
option.x0 = u0;
option.solver = 'CG';
[u,info] = mgHodgeLapE(Abar,f,node,elem,bdFlag,option); %#ok<*ASGLU>

Apmgoption.x0 = p0;
Apmgoption.printlevel = 1;
Apmgoption.freeDof = option.isFreeNode;
rg = g-G'*u;
[p,info2] = mg(Ap,rg,elem,Apmgoption);
% [p,info2] = amg(Ap,rg,Apmgoption);

u = u + grad*p;
p = -(DMinv*rg);

%% Output
time = cputime - t;
info.solverTime = time;

end