function pde = interfacedata(para)
if nargin == 0 
    betaMinus = 1;
    betaPlus = 10;
    c0 = -0.1;
    c1 = 0;
    omega = 6;
    xc = 0.02*sqrt(5);
    yc = 0.02*sqrt(5);
    r0 = 0.5;
    r1 = 0.2;
else
    if ~isstruct(para)
        exit('we need a struct data');
    end
    if ~isfield(para,'betaMinus') || isempty(para.betaMinus)
        betaMinus = 1;
    else
        betaMinus = para.betaMinus;
    end
    if ~isfield(para,'betaPlus') || isempty(para.betaPlus)
        betaPlus = 10;
    else
        betaPlus = para.betaPlus;
    end
    if ~isfield(para,'c0') || isempty(para.c0)
       c0 = -0.1;
    else
       c0 = para.c0;
    end
    if ~isfield(para,'c1') || isempty(para.c1)
        c1 = 0;
    else
        c1 = para.c1;
    end
    if ~isfield(para,'omega') || isempty(para.omega)
        omega = 12;
    else
        omega = para.omega;
    end
    if ~isfield(para,'xc') || isempty(para.xc)
        xc = 0.02*sqrt(5);
    else
        xc = para.xc;
    end
    if ~isfield(para,'yc') || isempty(para.yc)
        yc = 0.02*sqrt(5);
    else
        yc = para.yc;
    end
    if ~isfield(para,'r0') || isempty(para.r0)
        r0 = 0.5;
    else
        r0 = para.r0;
    end
    if ~isfield(para,'r1') || isempty(para.r1)
        r1 = 0.2;
    else
        r1 = para.r1;
    end
end

pde = struct('f',@f,'exactu',@exactu,'g_D',@g_D,'Du',@Du,'phi',@phi, 'exactw',@exactw,...
    'exactq',@exactq, 'd',@d,'exactuplus',@exactuplus,'exactuminus',@exactuminus,...
    'Duplus',@Duplus,'Duminus',@Duminus,'dplus',@dplus,'dminus',@dminus,'box',[-1,1,-1,1]);

% load data (right hand side function)
function z = d(p)
N = size(p,1);
z = zeros(N,1);
isInside = msign(phi(p)) <= 0;
z(isInside) = betaMinus;
z(~isInside) = betaPlus;
end
function z = dplus(p)
z = betaPlus;   
end
function z = dminus(p)
z = betaMinus;    
end
function z = f(p)
N = size(p,1);
z = zeros(N,1);
isInside = msign(phi(p))<=0;
rr = sum(p.^2,2);
z(isInside) = -4;
z(~isInside) = - 16*rr(~isInside);
end
% exact solution
function z = exactu(p)
N = size(p,1);
z = zeros(N,1);
isInside = msign(phi(p)) <= 0;
z(isInside) = exactuminus(p(isInside,:));
z(~isInside) = exactuplus(p(~isInside,:));
end

function z = exactuplus(p)
rr = sum(p.^2,2);
r = sqrt(rr);
rrrr = rr.^2;
r02 = r0*r0;
r04 = r02*r02;
z = (rrrr + c0 *log(2*r))/betaPlus + c1*(r02/betaMinus -(r04 + c0*log(2*r0))/betaPlus );
end
function z = exactuminus(p)
rr = sum(p.^2,2);
z = rr/betaMinus;
end

% Dirichlet boundary condition
function z = g_D(p)
z = exactuplus(p);
end
% Derivative of the exact solution
function z = Du(p)
N = size(p,1);
z = zeros(N,2);
isInside = msign(phi(p)) <= 0;
z(isInside,:) = Duminus(p(isInside,:));
z(~isInside,:) = Duplus(p(~isInside,:));
end

function z = Duminus(p)
z = 2*p/betaMinus;
end
function z = Duplus(p)
rr = sum(p.^2,2);
z(:,1) = (4*p(:,1).*rr + c0*p(:,1)./rr)/betaPlus;
z(:,2) = (4*p(:,2).*rr + c0*p(:,2)./rr)/betaPlus;
end

function z = exactw(p)
rr = sum(p.^2,2);
r = sqrt(rr);
rrrr = rr.^2;
r02 = r0*r0;
r04 = r02*r02;
z =  (rrrr + c0 *log(2*r))/betaPlus  +  c1*(r02/betaMinus -(r04 + c0*log(2*r0))/betaPlus );
z = z - rr/betaMinus;
end
function z = exactq(p)
% jump of the flux [\beta du/dn]


t = atan2((p(:,2) - yc),(p(:,1) - xc));
r = r0 + r1*sin(omega*t);
rp = omega*r1*cos(omega*t);
n1 = rp.*cos(t) - r.*sin(t);
n2 = rp.*sin(t) + r.*cos(t);
l = sqrt(n1.^2+n2.^2);
n = [n2./l, -n1./l];

rr = sum(p.^2,2);

z11 = (4*p(:,1).*rr + c0*p(:,1)./rr)/betaPlus;
z12 = (4*p(:,2).*rr + c0*p(:,2)./rr)/betaPlus;
z21 = 2*p(:,1)/betaMinus;
z22 = 2*p(:,2)/betaMinus;
z = [betaPlus*z11-betaMinus*z21, betaPlus*z12-betaMinus*z22];
z = sum(z.*n,2);
end

function z = phi(p) % level set function
theta = atan2((p(:,2) - yc),(p(:,1) - xc));
z = (p(:,2) - yc).^2 + (p(:,1)-xc).^2 - (r0 + r1*sin(omega*theta)).^2;
end

end