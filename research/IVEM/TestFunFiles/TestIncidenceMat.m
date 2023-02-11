%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form a face to edge incidence matrix

% faceDoFId = BasicFaceNum+(1:length(NewFaceID1));
% iief = reshape(repmat(faceDoFId,3,1),[],1);
% jjef = reshape(face2EDoF(NewFaceID1,:)',[],1);
% Sval = face2EDoFSign(NewFaceID1,:).*face2FDoFSign(NewFaceID1);
% xxef = reshape(Sval',[],1);
% ICface2edge = sparse(iief,jjef,xxef,NFdof,NEdof);
% CurlNew =  ICface2edge'*MatUI*ICface2edge;
% Dmatcurl = MatCurl - CurlNew;



BfeLocOriet1 = zeros(size(mesh.f_e(BasicFace,:),1),2);
BfeLocOriet2 = zeros(size(mesh.f_e(BasicFace,:),1),2);
BfeLocOriet3 = zeros(size(mesh.f_e(BasicFace,:),1),2);
BfN = mesh.f(BasicFace,:);
et1 = mesh.e(mesh.f_e(BasicFace,1),:);
et2 = mesh.e(mesh.f_e(BasicFace,2),:);
et3 = mesh.e(mesh.f_e(BasicFace,3),:);
for i = 1:3
    for j = 1:2
        IDj1 = et1(:,j); IDi1 = BfN(:,i);
        BfeLocOriet1(IDj1==IDi1,j) = i;
        IDj2 = et2(:,j); IDi2 = BfN(:,i);
        BfeLocOriet2(IDj2==IDi2,j) = i;
        IDj3 = et3(:,j); IDi3 = BfN(:,i);
        BfeLocOriet3(IDj3==IDi3,j) = i;
    end
end
BfeLocSign = zeros(size(mesh.f_e(BasicFace,:),1),3);
Pt1st = zeros(size(mesh.f_e(BasicFace,:),1),1);
for i =1:3
    index1 = min(abs(BfeLocOriet1 - i),[],2);
    index2 = min(abs(BfeLocOriet3 - i),[],2);
    IDi = (index1 == 0 & index2 == 0);
    Pt1st(IDi) = i;
end
ID1 = BfeLocOriet1(:,1) == Pt1st;
BfeLocSign(ID1,1) = 1; BfeLocSign(~ID1,1) = -1; 
ID3 = BfeLocOriet3(:,2) == Pt1st;
BfeLocSign(ID3,3) = 1; BfeLocSign(~ID3,3) = -1; 
ID2ndpt = zeros(size(mesh.f_e(BasicFace,:),1),1);
%%%%
Pt2nd = zeros(size(mesh.f_e(BasicFace,:),1),1);
index1 = BfeLocOriet1(:,1) - Pt1st;
Pt2nd(index1==0) = BfeLocOriet1(index1==0,2);
index2 = BfeLocOriet1(:,2) - Pt1st;
Pt2nd(index2==0) = BfeLocOriet1(index2==0,1);
ID2 = BfeLocOriet2(:,1) == Pt2nd;
BfeLocSign(ID2,2) = 1; BfeLocSign(~ID2,2) = -1;
%%%%%
Pt3rd = zeros(size(mesh.f_e(BasicFace,:),1),1);
index1 = BfeLocOriet2(:,1) - Pt2nd;
Pt3rd(index1==0) = BfeLocOriet2(index1==0,2);
index2 = BfeLocOriet2(:,2) - Pt2nd;
Pt3rd(index2==0) = BfeLocOriet2(index2==0,1);
%%%%
FaceNoriet = [Pt1st,Pt2nd,Pt3rd];
% f_e_orit = ones(size(mesh.f_e(BasicFace,:)));
% vtmp = TotalOldEdge(mesh.f_e(BasicFace,:));
% e_ind = [1,2; 2,3; 1,3];
% for i = 1:3
%     d_crect = gdof(vtmp(:,i),2) - gdof(vtmp(:,i),1);
%     f_e_orit(d_crect<0,i) = -1;
% end
% f_e_orit = f_e_orit.*Allface2FDoFSign(BasicFace);
et1 = node(Fgdof(1:BasicFaceNum,2),:) - node(Fgdof(1:BasicFaceNum,1),:);
et2 = node(Fgdof(1:BasicFaceNum,3),:) - node(Fgdof(1:BasicFaceNum,1),:);
AllFNormal = cross(et1,et2);
vOriet1 = [0,0,1]; vOriet2 = [1,0,0];
vid1 = sum(AllFNormal.*vOriet1,2)>10^(-12);
vid2 = (abs(sum(AllFNormal.*vOriet1,2))<=10^(-12) & sum(AllFNormal.*vOriet2,2)>10^(-12));
vid3 = (abs(sum(AllFNormal.*vOriet1,2))<=10^(-12) & abs(sum(AllFNormal.*vOriet2,2))<10^(-12)...
    & sum(AllFNormal.*vOriet3,2)>=10^(-12));
vid = (vid1+vid2+vid3>0);
Allface2FDoFSign = -ones(size(BfeLocOriet1,1),1);
Allface2FDoFSign(vid) = 1;
index1 = sum(abs(FaceNoriet - repmat([1,2,3],BasicFaceNum,1)),2)==0;
index2 = sum(abs(FaceNoriet - repmat([2,3,1],BasicFaceNum,1)),2)==0;
index3 = sum(abs(FaceNoriet - repmat([3,1,2],BasicFaceNum,1)),2)==0;
Signtemp = -ones(size(BfeLocOriet1,1),1);
index = index1+index2+index3>0;
Signtemp(index) = 1;
Allface2FDoFSign = Allface2FDoFSign.*Signtemp;

% faceDoFId = [1:BasicFaceNum,BasicFaceNum+NewFaceNintID'];
% iief = reshape(repmat(faceDoFId,3,1),[],1);
% vtmp = TotalOldEdge(mesh.f_e(BasicFace,:));
% jjef = reshape([vtmp;face2EDoF(NewFaceID1(NewFaceNintID),:)]',[],1);
% Sval = face2EDoFSign(NewFaceID1(NewFaceNintID),:).*face2FDoFSign(NewFaceID1(NewFaceNintID));
% xxef = reshape([BfeLocSign.*Allface2FDoFSign;Sval]',[],1);
% ICface2edge = sparse(iief,jjef,xxef,NFdof,NEdof);
% StiffCurlNew =  ICface2edge'*MassTotal*ICface2edge;
% Dmatcurl = StiffCurl - StiffCurlNew;

faceDoFId = 1:NFdof;
iief = reshape(repmat(faceDoFId,3,1),[],1);
vtmp = TotalOldEdge(mesh.f_e(BasicFace,:));
jjef = reshape([vtmp;face2EDoF(NewFaceID1,:)]',[],1);
Sval = face2EDoFSign(NewFaceID1,:).*face2FDoFSign(NewFaceID1);
% xxef = reshape([f_e_orit;Sval]',[],1);
xxef = reshape([BfeLocSign.*Allface2FDoFSign;Sval]',[],1);
ICface2edge = sparse(iief,jjef,xxef,NFdof,NEdof);
StiffCurlNew =  ICface2edge'*MassTotal*ICface2edge;
Dmatcurl = StiffCurl - StiffCurlNew;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form a local incidence matrix at one element
elemID = 1;
LocEDoF = currentElem(elemID,:); % for Hcurl
LocEDoFTMP = repmat(LocEDoF,size(Locelem2face2edge,1),1);
faceIDtmp = currentElem(1,:) - BasicFaceNum; % for Hdiv
Locelem2face2edge = face2EDoF(NewFaceID1(faceIDtmp),:);
LocIC = zeros(size(Locelem2face2edge,1),length(LocEDoF));
for i = 1:3
    LocICtmp = zeros(size(Locelem2face2edge,1),length(LocEDoF));
    for j = 1:size(LocEDoF,2)
        columi = Locelem2face2edge(:,i);
        columj = LocEDoFTMP(:,j);
        compareID = (columi==columj);
        LocICtmp(:,j) = compareID;
    end
    LocIC = LocIC + LocICtmp;
end
