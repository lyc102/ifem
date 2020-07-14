function [u,p,info,info2] =  mgMaxwellsaddle(A,G,f,g,node,elem,bdFlag,Me,grad,option,varargin)
%% mgMaxwellsaddle Solve the maxwell system with divgence free condition
%
%         [A  G] [u]  = [f]               
%         [G' O] [p]  = [g]               
%
% where  G = M_e*grad with the mass matrix for the edge element Me.
%
% This system can be rewritten as 
%         [A+G*DMinv*G'   G] [u]  = [f +G*DMinv*g]
%         [G'             O] [p]  = [g]
%
% where DMinv is the inverse of the diagonal matrix of the mass. Then
%
%  [A+G*DMinv*G'  G] [I    grad]                = [Abar   O  ]
%  [G'            O] [0  -Dminv*grad'*Me*grad ] = [G'     Ap]
%
% Then we compute the inverse Abar by mgHodgeLapE and Ap by mg. 
%
% Reference: 
% L Chen, Y Wu, L Zhong, J Zhou. MultiGrid Preconditioners for Mixed Finite
% Element Methods of the Vector Laplacian. Journal of Scientific Computing
% 77 (1), 101-128.
%
% Created by Long chen and Jie Zhou on Aug,2015.


Nf = length(f); 
Ng = length(g);
t = cputime;

%% Parameters
if ~exist('option','var'), option = []; end
Ndof = Nf + Ng; 
if ~isfield(option,'smoothingstep')  % smoothing steps
    option.smoothingstep = 3;
end
if ~isfield(option,'smoothingratio') % ratio of variable smoothing
    option.smoothingratio = 1.5;
end
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
if nargin>=6
    HB = varargin{1};
else
    HB = [];
end        
[u,info] = mgHodgeLapE(Abar,f,node,elem,bdFlag,option); %#ok<*ASGLU>
% [u,info] = mgHodgeLapE(Abar,f,node,elem,bdFlag,option); %#ok<*ASGLU>

Apoption.x0 = p0;
Apoption.printlevel = 1;
Apoption.freeDof = option.isFreeNode;
rg = g-G'*u;
[p,info2] = mg(Ap,rg,elem,Apoption,HB);
% [p,info2] = amg(Ap,rg,Apmgoption);

u = u + grad*p;
p = -(DMinv*rg);

%% Output
time = cputime - t;
info.solverTime = time;

end