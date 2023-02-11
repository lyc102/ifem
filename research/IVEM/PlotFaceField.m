Fid = 519;
Fnd = mesh.p(mesh.f(Fid,:),:);
idtmp = [1,2,3,1];
intpt0ID = mesh.eLoc(mesh.f_e(Fid,:));
intptID = -mesh.eLoc(mesh.f_e(Fid,intpt0ID<0));
intpt = mesh.eIntP(intptID,:);

plot3(Fnd(idtmp,1),Fnd(idtmp,2),Fnd(idtmp,3),'k','LineWidth',1)
hold on
plot3(intpt(:,1),intpt(:,2),intpt(:,3),'r','LineWidth',1)

testE1 = mesh.p(mesh.e(mesh.f_e(Fid,1),2),:) - intpt(1,:); % the first cutting edge corresponds to the first edge
mpt1 = 1/2*(mesh.p(mesh.e(mesh.f_e(Fid,1),2),:) + intpt(1,:));
testE2 = intpt(1,:) - mesh.p(mesh.e(mesh.f_e(Fid,1),1),:);
mpt2 = 1/2*(intpt(1,:) + mesh.p(mesh.e(mesh.f_e(Fid,1),1),:));
lE = norm(testE1) + norm(testE2);
ftId = -mesh.tLoc(mesh.f_t(Fid,:));
bas11 = squeeze(femI.bas1(ftId(1),:,:));
bas12 = squeeze(femI.bas2(ftId(1),:,:));
bas21 = squeeze(femI.bas1(ftId(2),:,:));
bas22 = squeeze(femI.bas2(ftId(2),:,:));

% use this plot to check the piece
% scatter3(mpt1(:,1),mpt1(:,2),mpt1(:,3),'r')

dof1 = zeros(6,1);
for i = 1:6
    
    dof1(i) = (dot((cross(bas11(1:3,i),mpt2') + bas11(4:6,i)),testE2)+ ...
        dot((cross(bas12(1:3,i),mpt1') + bas12(4:6,i)),testE1));
    
end

dof2 = zeros(6,1);
for i = 1:6
    
    dof2(i) = (dot((cross(bas21(1:3,i),mpt2') + bas21(4:6,i)),testE2)+ ...
        dot((cross(bas22(1:3,i),mpt1') + bas22(4:6,i)),testE1));
    
end
% can check which this edge is associated which dof, and which piece is
% related to bas1 or bas2
i1 = find(abs(dof1-1)<10^(-10)); i2 = find(abs(dof2-1)<10^(-10));
bas11 = bas11(:,i1); bas12 = bas12(:,i1);
bas21 = bas21(:,i2); bas22 = bas22(:,i2);

xm = min(Fnd(idtmp,1));
xM = max(Fnd(idtmp,1));
ym = min(Fnd(idtmp,2));
yM = max(Fnd(idtmp,2));
h = xM-xm;
[x,y] = meshgrid(xm:h/16:xM,ym:h/16:yM);
z = ones(size(x))*Fnd(1,3); % assume this face is parallel to xy plane
szh = size(x,1);
funx1 = zeros(size(x));
funy1 = zeros(size(x));
funz1 = zeros(size(x));

for i = 1:szh
    for j = 1:szh
        
        if x(i,j) < intpt(1,1)
            piece_ind = 1;
        elseif x(i,j) >= intpt(1,1)
            piece_ind = 2;
        end
        
        vec = ifebas(x(i,j), y(i,j),z(i,j),bas11,bas12,piece_ind);
        
        funx1(i,j) = vec(1);      
        funy1(i,j) = vec(2);     
        funz1(i,j) = vec(3);     
        
        if y(i,j)-x(i,j)>Fnd(1,2)-Fnd(1,1)
            funx1(i,j) = 0; funy1(i,j) = 0; funz1(i,j) = 0; 
        end
        
    end
end
quiver3(x,y,z+0.03,funx1,funy1,funz1,2,'b')


funx2 = zeros(size(x));
funy2 = zeros(size(x));
funz2 = zeros(size(x));

for i = 1:szh
    for j = 1:szh
        
        if x(i,j) < intpt(1,1)
            piece_ind = 1;
        elseif x(i,j) >= intpt(1,1)
            piece_ind = 2;
        end
        
        vec = ifebas(x(i,j), y(i,j),z(i,j),bas21,bas22,piece_ind);
        
        funx2(i,j) = vec(1);      
        funy2(i,j) = vec(2);     
        funz2(i,j) = vec(3);     
        
        if y(i,j)-x(i,j)>Fnd(1,2)-Fnd(1,1)
            funx2(i,j) = 0; funy2(i,j) = 0; funz2(i,j) = 0; 
        end
        
    end
end
quiver3(x,y,z-0.03,funx2,funy2,funz2,2,'color',[1,0.5,0])

box on
axis equal

Elem1Nd = mesh.p(mesh.t(mesh.f_t(Fid,1),:),:);
Elem2Nd = mesh.p(mesh.t(mesh.f_t(Fid,2),:),:);
idtmp = [1,2,3,1,4,3,1,4,2];
plot3(Elem1Nd(idtmp,1),Elem1Nd(idtmp,2),Elem1Nd(idtmp,3),'k','LineWidth',1)
hold on
plot3(Elem2Nd(idtmp,1),Elem2Nd(idtmp,2),Elem2Nd(idtmp,3),'k','LineWidth',1)
Eintpt1Id = mesh.eLoc(mesh.t_e(mesh.f_t(Fid,1),:));
Eintpt1Id = -mesh.eLoc(mesh.t_e(mesh.f_t(Fid,1),Eintpt1Id<0));
Eintpt1 = mesh.eIntP(Eintpt1Id,:);
Eintpt2Id = mesh.eLoc(mesh.t_e(mesh.f_t(Fid,2),:));
Eintpt2Id = -mesh.eLoc(mesh.t_e(mesh.f_t(Fid,2),Eintpt2Id<0));
Eintpt2 = mesh.eIntP(Eintpt2Id,:);
patch(Eintpt1(:,1),Eintpt1(:,2),Eintpt1(:,3),'r','FaceAlpha',0.7)
idtmp = [1,2,4,3];
patch(Eintpt2(idtmp,1),Eintpt2(idtmp,2),Eintpt2(idtmp,3),'r','FaceAlpha',0.7)
box on
axis off
axis equal


function u = ifebas(x,y,z,bas1,bas2,piecie)
fnorm = [0;0;1]; % assume this face is parallel to xy plane

if piecie == 1
    bas = bas1;
elseif piecie == 2
    bas = bas2;
end

u = cross(bas(1:3),[x;y;z]) + bas(4:6);
u = cross(u,fnorm);

end

