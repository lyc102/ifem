function e = mgVcyclefracLap2d(A,r,pde,cube,h,s,option)

global mA

%% Setting
% smoother parameters
if isfield(option,'smootherType')
    smootherType = option.smootherType;
else
    smootherType = 3; 
end
if smootherType == 2 % Jacobi smoother
    if isfield(option,'omega')
        omega = option.omega;
    else
        omega = 0.5;
    end
end
if isfield(option,'smoothingstep')
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
y0 = cube(5); y1 = cube(6);
nx1 = (cube(2)-cube(1))/h+1;  % number of grid points in x1-direction
nx2 = (cube(4)-cube(3))/h+1;  % number of grid points in x2-direction
ny = (y1-y0)/h+1;  % number of grid points in y-direction
nx = nx1*nx2;
% nynx1 = ny*nx1;
% 1D vector to 3D matrix
r3d = reshape(r,ny,nx1,nx2);
e = zeros(N,1);
e3d = reshape(e,ny,nx1,nx2);
% index for Red-Black smoothing
nodeidx = reshape(1:N,ny,nx1,nx2);
nodeidx2d = reshape(nodeidx,ny,nx);
c1idx = nodeidx2d(2:2:nx1-1,2:2:nx2-1); 
c2idx = nodeidx2d(3:2:nx1-2,2:2:nx2-1);
c3idx = nodeidx2d(2:2:nx1-1,3:2:nx2-2);
c4idx = nodeidx2d(3:2:nx1-2,3:2:nx2-2);
c1linear = nodeidx2d(:,c1idx);
c2linear = nodeidx2d(:,c2idx);
c3linear = nodeidx2d(:,c3idx);
c4linear = nodeidx2d(:,c4idx);
for it = 1:smoothingstep
    switch smootherType
        case 0  % Point G-S Smoothing
            e = e + tril(A)\r;
            r = rold - A*e;
            e3d = reshape(e,ny,nx1,nx2);
        case 1
            for k1 = 2:nx1-1 % only update interiori lines
                for k2 = 2:nx2-1 % only update interiori lines
                    % form residual
%                     idx = (k2-1)*nynx1+(k1-1)*ny+(1:ny-1);
                    idx = nodeidx(:,k1,k2);
                    de = A(idx,idx)\r3d(:,k1,k2);
                    e3d(:,k1,k2) = e3d(:,k1,k2) + de;
                    r3d(:,k1,k2) = 0;  % exact solve
                    % update neighbor columns
                    r3d(:,k1+1,k2) = r3d(:,k1+1,k2) - A(nodeidx(:,k1+1,k2),idx)*de;
                    r3d(:,k1-1,k2) = r3d(:,k1-1,k2) - A(nodeidx(:,k1-1,k2),idx)*de;
                    r3d(:,k1,k2+1) = r3d(:,k1,k2+1) - A(nodeidx(:,k1,k2+1),idx)*de;
                    r3d(:,k1,k2-1) = r3d(:,k1,k2-1) - A(nodeidx(:,k1,k2-1),idx)*de;
                    r3d(:,k1+1,k2+1) = r3d(:,k1+1,k2+1) - A(nodeidx(:,k1+1,k2+1),idx)*de;
                    r3d(:,k1-1,k2-1) = r3d(:,k1-1,k2-1) - A(nodeidx(:,k1-1,k2-1),idx)*de;
                    r3d(:,k1-1,k2+1) = r3d(:,k1-1,k2+1) - A(nodeidx(:,k1-1,k2+1),idx)*de;
                    r3d(:,k1+1,k2-1) = r3d(:,k1+1,k2-1) - A(nodeidx(:,k1+1,k2-1),idx)*de;
                end
            end
        case 2 % weighted Jacobi
            idx = nodeidx(:,3,3); % take one typical column            
            r2d = reshape(r,ny,nx1*nx2); % 3-D matrix to a 2-D matrix            
            de = A(idx,idx)\r2d;
            e3d = e3d + omega*reshape(de,ny,nx1,nx2);
            r = rold - A*(e3d(:));
            r3d = reshape(r,ny,nx1,nx2);
        case 3 % four colors Gauss-Seidel
            idx = nodeidx(:,2,2); % take one typical column
%             r2d = reshape(r,ny,nx1*nx2); % 3-D matrix to a 2-D matrix
            r2d = reshape(r3d,ny,nx1*nx2);
            e2d = reshape(e3d,ny,nx1*nx2);
            % color 1
            de = A(idx,idx)\r2d(:,c1idx);
            e2d(:,c1idx) = e2d(:,c1idx) + de;
            % update residual
            r2d(:,c1idx) = 0;
            neigidx = [c2idx(:); c3idx(:); c4idx(:)];
            neigidxlin = [c2linear(:); c3linear(:); c4linear(:)];
            r2d(:,neigidx) = r2d(:,neigidx)-reshape(A(neigidxlin,c1linear)*de(:),ny,length(neigidx));
            % previous r2d(:,blackidx) = 0;
            % color 2
            de = A(idx,idx)\r2d(:,c2idx);
            e2d(:,c2idx) = e2d(:,c2idx) + de;
            % update residual
            r2d(:,c2idx) = 0;
            neigidx = [c1idx(:); c3idx(:); c4idx(:)];
            neigidxlin = [c1linear(:); c3linear(:); c4linear(:)];
            r2d(:,neigidx) = r2d(:,neigidx)-reshape(A(neigidxlin,c2linear)*de(:),ny,length(neigidx));
            % color 3
            de = A(idx,idx)\r2d(:,c3idx);
            e2d(:,c3idx) = e2d(:,c3idx) + de;
            % update residual
            r2d(:,c3idx) = 0;
            neigidx = [c1idx(:); c2idx(:); c4idx(:)];
            neigidxlin = [c1linear(:); c2linear(:); c4linear(:)];
            r2d(:,neigidx) = r2d(:,neigidx)-reshape(A(neigidxlin,c3linear)*de(:),ny,length(neigidx));
            % color 4
            de = A(idx,idx)\r2d(:,c4idx);
            e2d(:,c4idx) = e2d(:,c4idx) + de;
            % update residual
            r2d(:,c4idx) = 0;
            neigidx = [c1idx(:); c2idx(:); c3idx(:)];
            neigidxlin = [c1linear(:); c2linear(:); c3linear(:)];
            r2d(:,neigidx) = r2d(:,neigidx)-reshape(A(neigidxlin,c4linear)*de(:),ny,length(neigidx));
            % previous r2d(:,redidx) = 0;            
            % change back to single index
            r3d = reshape(r2d,ny,nx1,nx2);
            e3d = reshape(e2d,ny,nx1,nx2);
    end
end

%% Transfer operator 
% prolongation and restriction in x-direction
Ix1 = prolongation1d(log2(nx1-1)-1);
Rx1 = Ix1';
Ix2 = prolongation1d(log2(nx2-1)-1);
Rx2 = Ix2';
% prolongation and restriction in y-direction
nx1c = (nx1 - 1)/2 + 1;
nx2c = (nx2 - 1)/2 + 1;
% geometric quantity
nyc = (ny - 1)/2 + 1;
% My = ny - 1; 
Tyf = gradmap(y0,y1,gamma,h);
% Tyf = ((0:My)'/My).^gamma*(y1-y0) + y0;
hf = diff(Tyf);
% Myc = nyc - 1;
Tyc = gradmap(y0,y1,gamma,2*h);
% Tyc = ((0:Myc)'/Myc).^gamma*(y1-y0) + y0;
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
rc3d = zeros(nyc,nx1c,nx2c);
% rc3d = tprod(Ry,[1,-1],r3d,[-1,2,3]); % restriction in y-direction
% rc3d = tprod(Rx1,[2,-1],rc3d,[1,-1,3]);
% rc3d = tprod(Rx2,[3,-1],rc3d,[1,2,-1]);
rc3dtemp = reshape(Ry*reshape(r3d,ny,nx1*nx2),nyc,nx1,nx2); % restriction in y-direction
for iy = 1:nyc
    rc3d(iy,:,:) = Rx1*squeeze(rc3dtemp(iy,:,:))*Ix2;
end
rc = rc3d(:);

%% Coarse grid correction
level = -log2(h);
Ac = mA{level-1};
if level <= 3
    ec = Ac\rc;
else
    ec = mgVcyclefracLap2d(Ac,rc,pde,cube,2*h,s,option);
end

%% Prolongation
ec3d = reshape(ec,nyc,nx1c,nx2c);
ec3dtemp = zeros(nyc,nx1,nx2);
for iy = 1:nyc
    ec3dtemp(iy,:,:) = Ix1*squeeze(ec3d(iy,:,:))*Rx2;
end
ec3d = reshape(Iy*reshape(ec3dtemp,nyc,nx1*nx2),ny,nx1,nx2); % prolongation in y-direction
% ec3d = tprod(Iy,[1,-1],ec3d,[-1,2,3]);
% ec3d = tprod(Ix1,[2,-1],ec3d,[1,-1,3]);
% ec3d = tprod(Ix2,[3,-1],ec3d,[1,2,-1]);
e3d = e3d + ec3d;
e = e3d(:);
% update residual when e updated. remember we are solving Ae = rold.
r = rold - A*(e3d(:));
r3d = reshape(r,ny,nx1,nx2);

%% Post-smoothing
for it = 1:smoothingstep
    switch smootherType
        case 0 % pointwise G-S smoothing
            e = e + triu(A)\r;
            r = rold - A*e;
            e3d = reshape(e,ny,nx1,nx2);
        case 1
            for k1 = nx1-1:-1:2 % only update interiori lines
                for k2 = nx2-1:-1:2 % only update interiori lines
                    % form residual
                    idx = nodeidx(:,k1,k2);                    
                    de = A(idx,idx)\r3d(:,k1,k2);
                    e3d(:,k1,k2) = e3d(:,k1,k2) + de;
                    r3d(:,k1,k2) = 0;  % exact solve
                    % update neighbor columns
                    r3d(:,k1+1,k2) = r3d(:,k1+1,k2) - A(nodeidx(:,k1+1,k2),idx)*de;
                    r3d(:,k1-1,k2) = r3d(:,k1-1,k2) - A(nodeidx(:,k1-1,k2),idx)*de;
                    r3d(:,k1,k2+1) = r3d(:,k1,k2+1) - A(nodeidx(:,k1,k2+1),idx)*de;
                    r3d(:,k1,k2-1) = r3d(:,k1,k2-1) - A(nodeidx(:,k1,k2-1),idx)*de;
                    r3d(:,k1+1,k2+1) = r3d(:,k1+1,k2+1) - A(nodeidx(:,k1+1,k2+1),idx)*de;
                    r3d(:,k1-1,k2-1) = r3d(:,k1-1,k2-1) - A(nodeidx(:,k1-1,k2-1),idx)*de;
                    r3d(:,k1-1,k2+1) = r3d(:,k1-1,k2+1) - A(nodeidx(:,k1-1,k2+1),idx)*de;
                    r3d(:,k1+1,k2-1) = r3d(:,k1+1,k2-1) - A(nodeidx(:,k1+1,k2-1),idx)*de;
                end
            end
        case 2 % weighted Jacobi
            idx = nodeidx(:,3,3); % take one typical column            
            r2d = reshape(r,ny,nx1*nx2); % 3-D matrix to a 2-D matrix            
            de = A(idx,idx)\r2d;
            e3d = e3d + omega*reshape(de,ny,nx1,nx2);
            r = rold - A*(e3d(:));
            r3d = reshape(r,ny,nx1,nx2);
        case 3 % red-black block Gauss-Seidel
            idx = nodeidx(:,2,2); % take one typical column
%             r2d = reshape(r,ny,nx1*nx2); % 3-D matrix to a 2-D matrix
            r2d = reshape(r3d,ny,nx1*nx2);
            e2d = reshape(e3d,ny,nx1*nx2);
            % color 1
            de = A(idx,idx)\r2d(:,c1idx);
            e2d(:,c1idx) = e2d(:,c1idx) + de;
            % update residual
            r2d(:,c1idx) = 0;
            neigidx = [c2idx(:); c3idx(:); c4idx(:)];
            neigidxlin = [c2linear(:); c3linear(:); c4linear(:)];
            r2d(:,neigidx) = r2d(:,neigidx)-reshape(A(neigidxlin,c1linear)*de(:),ny,length(neigidx));
            % previous r2d(:,blackidx) = 0;
            % color 2
            de = A(idx,idx)\r2d(:,c2idx);
            e2d(:,c2idx) = e2d(:,c2idx) + de;
            % update residual
            r2d(:,c2idx) = 0;
            neigidx = [c1idx(:); c3idx(:); c4idx(:)];
            neigidxlin = [c1linear(:); c3linear(:); c4linear(:)];
            r2d(:,neigidx) = r2d(:,neigidx)-reshape(A(neigidxlin,c2linear)*de(:),ny,length(neigidx));
            % color 3
            de = A(idx,idx)\r2d(:,c3idx);
            e2d(:,c3idx) = e2d(:,c3idx) + de;
            % update residual
            r2d(:,c3idx) = 0;
            neigidx = [c1idx(:); c2idx(:); c4idx(:)];
            neigidxlin = [c1linear(:); c2linear(:); c4linear(:)];
            r2d(:,neigidx) = r2d(:,neigidx)-reshape(A(neigidxlin,c3linear)*de(:),ny,length(neigidx));
            % color 4
            de = A(idx,idx)\r2d(:,c4idx);
            e2d(:,c4idx) = e2d(:,c4idx) + de;
            % update residual
            r2d(:,c4idx) = 0;
            neigidx = [c1idx(:); c2idx(:); c3idx(:)];
            neigidxlin = [c1linear(:); c2linear(:); c3linear(:)];
            r2d(:,neigidx) = r2d(:,neigidx)-reshape(A(neigidxlin,c4linear)*de(:),ny,length(neigidx));
            % previous r2d(:,redidx) = 0;            
            % change back to single index
            r3d = reshape(r2d,ny,nx1,nx2);
            e3d = reshape(e2d,ny,nx1,nx2);
    end
end
e = e3d(:);