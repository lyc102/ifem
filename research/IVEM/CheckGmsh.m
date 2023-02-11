BoxTorusMesh

% f = @(p) (sqrt(p(:,1).^2+p(:,2).^2)-0.6283).^2 + (p(:,3)+0.3).^2 - 0.2^2;
% mesh.t=Tet(:,1:4);
% mesh.p=Node(:,1:3);
% mesh.tLoc = Tet(:,5);
% tm = (mesh.p(mesh.t(:,1),:) + mesh.p(mesh.t(:,2),:) + mesh.p(mesh.t(:,3),:) + mesh.p(mesh.t(:,4),:))/4;
% fmval = feval(f,[tm(:,1),tm(:,2),tm(:,3)]);
% piece1IDm = find(fmval<0);
% piece2IDm = find(fmval>0);
% piece1ID = find(Tet(:,end)==2);
% piece2ID = find(Tet(:,end)==1);
% tet1 = mesh.t(piece1ID,:);
% T1 = auxstructure3(tet1);
% trisurf(T1.bdFace,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3))

mesh.t = Tetrahedra(:,1:4);
mesh.p = Vertices(:,1:3);
mesh.tLoc = 3-Tetrahedra(:,5);
tID1 = (mesh.tLoc == 1);
tID2 = (mesh.tLoc == 2);
T1 =  auxstructure3(mesh.t(tID1,:));
T2 =  auxstructure3(mesh.t(tID2,:));
NodeIntf2Idtmp = unique(reshape(T2.bdFace,[],1));
NodeIntf2tmp = mesh.p(NodeIntf2Idtmp,:);
idtmp = (NodeIntf2tmp(:,1)>-1+0.05).*(NodeIntf2tmp(:,1)<1-0.05).*...
    (NodeIntf2tmp(:,2)>-1+0.05).*(NodeIntf2tmp(:,2)<1-0.05).*...
    (NodeIntf2tmp(:,3)>-1+0.05).*(NodeIntf2tmp(:,3)<1-0.05);

NodeIntf1Id = unique(reshape(T1.bdFace,[],1));