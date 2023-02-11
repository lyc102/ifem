function pLoc = FindElemId(p,mesh)
% find which element that the point p belongs to
% p is a n*3 matrix of which each row is a point

x_min=min(mesh.p(:,1)); x_max=max(mesh.p(:,1));
y_min=min(mesh.p(:,2)); y_max=max(mesh.p(:,2));
z_min=min(mesh.p(:,3)); z_max=max(mesh.p(:,3));
hx = max(mesh.p(mesh.t(1,:),1)) - min(mesh.p(mesh.t(1,:),1));
hy = max(mesh.p(mesh.t(1,:),2)) - min(mesh.p(mesh.t(1,:),2));
hz = max(mesh.p(mesh.t(1,:),3)) - min(mesh.p(mesh.t(1,:),3));
nx = round((x_max-x_min)/hx);
ny = round((y_max-y_min)/hy);
nz = round((z_max-z_min)/hz);
n = nx*ny*nz;

x = p(:, 1) - x_min;
y = p(:, 2) - y_min;
z = p(:, 3) - z_min;

i = floor(x/hx)+1;
j = floor(y/hy)+1;
k = floor(z/hz)+1;

LocId = (k-1)*nx*ny + (j-1)*nx + i;
vert = mesh.T(LocId,:);
pref = p - [(i-1)*hx+x_min,(j-1)*hy+y_min,(k-1)*hz+z_min];
p0 = [0,0,0; hx,0,0; hx,hy,0; 0,hy,0;...
    0,0,hz; hx,0,hz; hx,hy,hz; 0,hy,hz];

tr = [1,2,3,7;
    1,6,2,7;
    1,5,6,7;
    1,8,5,7;
    1,4,8,7;
    1,3,4,7];

%pLocTmp = zeros(size(p,1),6);

for r = 1:6
    trjcoord = p0(tr(r,:),:);
    TR = triangulation(tr,p0);
    ID = pointLocation(TR,pref);
    %m = sum(trjcoord,1)/4;
    %trjcoordNew = (trjcoord - ones(4,1)*m)*10^(-7)+trjcoord;
    %ptInTriOrNot = funPtInTriCheckNd(trjcoord,pref);
%     if ptInTriOrNot == 1
%         pLoc(k) = (j-1)*n + LocInd(k);
%         
%         vk = mesh.p(mesh.t(pLoc(k),:),:);
%         uk = uh(mesh.t(pLoc(k),:));
%         T = vk(1:3,:) - ones(3,1)*vk(4,:); T = T';
%         bk = pk' - vk(4,:)';
%         lambda = T\bk;
%         uval(k) = lambda(1)*uk(1) + lambda(2)*uk(2) +...
%             lambda(3)*uk(3) + (1-sum(lambda))*uk(4);
%         
%     end
end

if min(ID)==0
    stp=1;
end
pLoc = (ID-1)*n + LocId;
    


%uval = zeros(size(p,1),1);

