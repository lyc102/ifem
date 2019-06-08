function e = mgVcyclefracLap1d(A,r,pde,square,h,s,option)

global mA

%% Setting
if ~isfield('option','smootherType')
    smootherType = option.smootherType;
else
    smootherType = 3; 
end
if ~isfield('option','smootherstep')
    smoothingstep = option.smoothingstep;
else
    smoothingstep = 3; 
end
if ~exist('smootherType','var'), smootherType = 3; end
% mapping parameter
if s == 0.5
    gamma = 1;
else
    gamma = 3/(2*s)+0.1;
end

%% line smoothing
N = length(r);
rold = r;
x0 = square(1); x1 = square(2); 
y0 = square(3); y1 = square(4);
nx = (x1-x0)/h+1;  % number of grid points in x-direction
ny = (y1-y0)/h+1;  % number of grid points in x-direction
% 1D vector to 2D matrix
r2d = reshape(r,ny,nx);
e = zeros(N,1);
e2d = reshape(e,ny,nx);
% index for Red-Black smoothing
nodeidx = reshape(1:N,ny,nx);
redidx = 2:2:nx-1;
blackidx = 3:2:nx-2; 
rlinear = nodeidx(:,redidx);
blinear = nodeidx(:,blackidx);
for it = 1:smoothingstep
    switch smootherType
        case 0
            e = e + tril(A)\r;
            r = rold - A*e;
            e2d = reshape(e,ny,nx);            
        case 1
            for k = 2:nx-1 % only update interiori lines
                % form residual
                idx = (k-1)*ny+1:k*ny;  % the end nodes is included?
                de = A(idx,idx)\r2d(:,k);
                e2d(:,k) = e2d(:,k) + de;
                r2d(:,k) = 0;  % exact solve
                if k<nx-1   % update right column
                    r2d(:,k+1) = r2d(:,k+1) - A(idx+ny,idx)*de;
                end
                if k>2    % update left column
                % for left-to-right ordering, the left r(idx-ny) is zero
                    r2d(:,k-1) = - A(idx-ny,idx)*de;
                end
            end
        case 2 % weighted Jacobi
            for k = 2:nx-1
                % form residual
                idx = (k-1)*ny+1:k*ny;
                de = A(idx,idx)\r(idx);
                e(idx) = e(idx) + 0.25*de;
            end            
            r = rold - A*e;
        case 3 % red-black block Gauss-Seidel
            idx = ny+1:2*ny;
            % 1: red lines
            de = A(idx,idx)\r2d(:,redidx);
            e2d(:,redidx) = e2d(:,redidx) + de;
            % update residual
            r2d(:,redidx) = 0;
            r2d(:,blackidx) = r2d(:,blackidx)-reshape(A(blinear,rlinear)*de(:),ny,length(blackidx));
            % previous r2d(:,blackidx) = 0;
            % 2: black lines
            de = A(idx,idx)\r2d(:,blackidx);
            e2d(:,blackidx) = e2d(:,blackidx) + de;
            % update residual
            r2d(:,blackidx) = 0;
            r2d(:,redidx) = -reshape(A(rlinear,blinear)*de(:),ny,length(redidx));
            % previous r2d(:,redidx) = 0;            
    end
end

%% Transfer operator 
% prolongation and restriction in x-direction
Ix = prolongation1d(log2(nx-1)-1);
Rx = Ix';
% prolongation and restriction in y-direction
nxc = (nx - 1)/2 + 1;
nyc = (ny - 1)/2 + 1;
% geometric quantity
My = ny - 1; 
Tyf = ((0:My)'/My).^gamma*(y1-y0) + y0;
hf = diff(Tyf);
Myc = nyc - 1;
Tyc = ((0:Myc)'/Myc).^gamma*(y1-y0) + y0;
hc = diff(Tyc);
jc = 2:nyc-1;  % interiori points in the coarse grid
j = 2*jc-1;  % index of coarse points in the fine grids
alpha = zeros(1,nyc);
beta = zeros(1,nyc);
alpha(jc) = hf(j)./hc(jc);
beta(jc) = hf(j-1)./hc(jc-1);
alpha(1) = 1-beta(2);
jc = 2:(nyc-1);
jf = 2*jc-1;  % index of coarse points in the fine grids
% due to the Neumann boundary condition, 1 is included
ii = [1 jf 2 jf+1 jf-1];
jj = [1 jc 1 jc jc];
ss = [ones(1,nyc-1) alpha(1:nyc-1) beta(2:nyc-1)];
Iy = sparse(ii,jj,ss,ny,nyc);
Ry = Iy';

%% Restriction
% 1D vector to 2D matrix
% r2d = reshape(r,ny,nx);
rc2d = Ry*r2d*Ix;
rc = rc2d(:);

%% Coarse grid correction
% option.solver = 'none';
% [u,eqn] = fracLap1d(square,2*h,pde,option);
level = -log2(h);
Ac = mA{level-1};
fixedNode = [1:nyc (2:nxc-1)*nyc (nxc-1)*nyc+(1:nyc)];
isBdNode = false(length(rc),1);
isBdNode(fixedNode) = true;
freeNode = find(~isBdNode);
if level <= 3
    ec = zeros(length(rc),1);
    ec(freeNode) = Ac(freeNode,freeNode)\rc(freeNode);
else
    ec = mgVcyclefracLap1d(Ac,rc,pde,square,2*h,s,option);
end

%% Prolongation
% 1D vector to 2D matrix
ec2d = reshape(ec,nyc,nxc);
e2d = e2d + Iy*ec2d*Rx;
e = e2d(:);
% update residual when e updated. remember we are solving Ae = rold.
r = rold - A*(e2d(:));
r2d = reshape(r,ny,nx);

%% Post-smoothing
for it = 1:smoothingstep
    switch smootherType
        case 0 % pointwise G-S smoothing
            e = e + triu(A)\r;
            r = rold - A*e;
            e2d = reshape(e,ny,nx);
        case 1
            for k = nx-1:-1:2
                % form residual
                idx = (k-1)*ny+1:k*ny;
                de = A(idx,idx)\r2d(:,k);
                e2d(:,k) = e2d(:,k) + de;
                r2d(:,k) = 0;                
                if k < nx-1   % update right column
                % for right-to-left ordering, the right r(idx+ny) is zero
                    r2d(:,k+1) = - A(idx+ny,idx)*de;
                end
                if k>2    % update previous column
                    r2d(:,k-1) = r2d(:,k-1) - A(idx-ny,idx)*de;
                end
            end   
        case 2 % weighted Jacobi
            for k = 1:nx
                % form residual
                idx = (k-1)*ny+1:k*ny;
                de = A(idx,idx)\r(idx);
                e(idx) = e(idx) + 0.25*de;
            end            
            r = rold - A*e;            
        case 3 % red-black block Gauss-Seidel
            idx = ny+1:2*ny;
            % 1: red lines
            de = A(idx,idx)\r2d(:,redidx);
            e2d(:,redidx) = e2d(:,redidx) + de;
            % update residual
            r2d(:,redidx) = 0;
            r2d(:,blackidx) = r2d(:,blackidx)-reshape(A(blinear,rlinear)*de(:),ny,length(blackidx));
            % previous r2d(:,blackidx) = 0;
            % 2: black lines
            de = A(idx,idx)\r2d(:,blackidx);
            e2d(:,blackidx) = e2d(:,blackidx) + de;
            % update residual
            r2d(:,blackidx) = 0;
            r2d(:,redidx) = -reshape(A(rlinear,blinear)*de(:),ny,length(redidx));
            % previous r2d(:,redidx) = 0;            
    end
end
e = e2d(:);