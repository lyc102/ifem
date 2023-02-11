function mesh = genIntfMesh3D(mesh,intf)

%% USAGE: enrich mesh information for interface problem.
% INPUTS: 
% mesh --- nt-by-3 matrix contains p and t information
% intf --- interface function in analytical form
%               g(x,y,z) = 0
%
% OUTPUTS:
% pLoc --- mesh node indices.  possible values:
%               pLoc(i) =-1 if g(xi,yi,zi) < 0 
%               pLoc(i) = 0 if g(xi,yi,zi) = 0 (on interface) 
%               pLoc(i) = 1 if g(xi,yi,zi) > 0 
% tLoc --- mesh element indices.  possible values:
%               tLoc(i) = 1 if the element Ti is in region 1 (g(X) < 0)
%               tLoc(i) = 2 if the element Ti is in region 2 (g(X) > 0)
%               tLoc(i) = -ki if Ti is the ki-th interface element
% intType --- interface element type.  possible values:
%             0 for noninterface element, 
%             1 for Type 1: one vertex on one side, three vertices on the other
%             2 for Type 2: two vertice on each side of interface)
% eLoc --- mesh edge indices.  possible values:
%               eLoc(i) = 1 if the element Ti is in region 1 (g(X) < 0)
%               eLoc(i) = 2 if the element Ti is in region 2 (g(X) > 0)
%               eLoc(i) = -ki if Ti is the ki-th interface element
% eIntP --- neI-by-3 vector, stores the coordinate of interface
%                    intersection points on each edge. 
%           
% Last updated by Xu Zhang on 07/09/20

%% 1. generate pLoc
%% basic assumption: 
p = mesh.p; t = mesh.t; e = mesh.e; f = mesh.f;
nt = length(t); ne = length(e); nf = length(f);
pLoc = sign(feval(intf,p(:,1),p(:,2),p(:,3))); pLoc(pLoc == 0) = 1;

%% 2. generate tLoc and intType
id1 = sign(feval(intf,p(t(:,1),1),p(t(:,1),2),p(t(:,1),3)));
id2 = sign(feval(intf,p(t(:,2),1),p(t(:,2),2),p(t(:,2),3)));
id3 = sign(feval(intf,p(t(:,3),1),p(t(:,3),2),p(t(:,3),3)));
id4 = sign(feval(intf,p(t(:,4),1),p(t(:,4),2),p(t(:,4),3)));
id = id1+id2+id3+id4;
idm = id1.*id2.*id3.*id4;
tLoc = zeros(nt,1); intType = zeros(nt,1);
r1 = find(id == -4); tLoc(r1) = 1; intType(r1) = 0;
r2 = find(id == 4);  tLoc(r2) = 2; intType(r2) = 0;
r6 = (id == -1); tLoc(r6) = 1; intType(r6) = 0; % one pt is - and others are zero
r7 = (id == 1); tLoc(r7) = 2;  intType(r7) = 0; % one pt is + and others are zero
r8 = (id == -2 & idm == 0); tLoc(r8) = 1; intType(r8) = 0; % two pts are - and others are zero
r9 = (id == 2 & idm == 0);  tLoc(r9) = 2; intType(r8) = 0; % two pts are + and others are zero
r10 = (id == -3 & idm == 0); tLoc(r10) = 1; intType(r10) = 0; % three pts are - and the other s zero
r11 = (id == 3 & idm == 0);  tLoc(r11) = 2; intType(r11) = 0; % three pts are + and the other s zero
r3 = (id == -2 & idm ~= 0); intType(r3) = 1;
r4 = (id == 2 & idm ~= 0);  intType(r4) = 1;
r5 = (id == 0); intType(r5) = 2;
nt1 = length(r1) + sum(r6) + sum(r8) + sum(r10); 
nt2 = length(r2) + sum(r7) + sum(r9) + sum(r11); 
ntI = nt-nt1-nt2;
tLoc(tLoc == 0) = -(1:ntI);

%% 3. generate eLoc
eLoc = zeros(ne,1); 
ide1 = sign(feval(intf,p(e(:,1),1),p(e(:,1),2),p(e(:,1),3))); ide1(ide1 == 0) = 1;
ide2 = sign(feval(intf,p(e(:,2),1),p(e(:,2),2),p(e(:,2),3))); ide2(ide2 == 0) = 1;
ide = ide1+ide2;
re1 = find(ide == -2); eLoc(re1) = 1; 
re2 = find(ide == 2);  eLoc(re2) = 2; 
ne1 = length(re1); ne2 = length(re2); neI = ne-ne1-ne2;
eLoc(eLoc == 0) = -(1:neI);

%% 4. generate eIntP
p1 = p(e(eLoc<0,1),:);
p2 = p(e(eLoc<0,2),:);
eIntP = IntersectPoint3D(p1, p2, intf);

%% 5. generate fLoc
fLoc = zeros(nf,1); 
fID = sum(pLoc(mesh.f),2);
rf1 = find(fID == -3); fLoc(rf1) = 1;
rf2 = find(fID == 3);  fLoc(rf2) = 2;
nf1 = length(rf1); nf2 = length(rf2); nfI = nf-nf1-nf2;
fLoc(fLoc == 0) = -(1:nfI);

%% Form Mesh
mesh.pLoc = pLoc;
mesh.eLoc = eLoc;
mesh.tLoc = tLoc;
mesh.fLoc = fLoc;
mesh.eIntP = eIntP;
mesh.intType = intType;
