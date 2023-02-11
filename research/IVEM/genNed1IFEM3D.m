function femI = genNed1IFEM3D(mesh,fem,btM,btP,ahM,ahP,option)
%% Usage: Generate Quadrature Information on Interface Elements
%         Each interface tetrahedron is cut into 4 to 9 small tetrahedra 
%
% femI.t --- global index on small tetrahedra
% femI.bas --- nt-6-6 matrix stores basis functions on all small tetra   
% femI.gx --- Gaussian x nodes
% femI.gy --- Gaussian y nodes
% femI.gz --- Gaussian z nodes
% femI.gw --- Gaussian weights
% femI.area --- Areas of all small tetra
% femI.quadT --- number of small tetras of each interface element
% femI.locID --- order of local index, to match with the basis function

% ***** The following info only used for partial penalty IFEM *****

% femI.bas1 --- nti-6-6 basis function on 1st piece of interface element
% femI.bas2 --- nti-6-6 basis function on 2nd piece of interface element
% femI.plusPC --- 1 or 2 denotes which piece corresponds to btP.

% Last Modified by Ruchi Guo 19/11/2020

%% 
ntI = -min(mesh.tLoc);
g2ldof = zeros(8*ntI,6); bas = zeros(8*ntI,6,6); A = zeros(8*ntI,1); tIntfID = zeros(8*ntI,1);
if nargin == 6
    %bcind = [1,1,1,1,1,1];  
    ng = 4;
end
if nargin == 7
    ng = option.ng;
end
gx = zeros(8*ntI,ng); gy = zeros(8*ntI,ng); gz = zeros(8*ntI,ng);
quadT = zeros(ntI,1); locIDvec = zeros(ntI,4); elocIDvec = zeros(ntI,6);
bas1 = zeros(ntI,6,6); bas2 = zeros(ntI,6,6); plusPC = zeros(ntI,1);
intID = find(mesh.tLoc<0);
tetID = zeros(ntI,2);
id = 0;


for i = 1:ntI
    tID = intID(i);
    t_e = mesh.t_e(tID,:);
    intpt0 = mesh.eIntP(-mesh.eLoc(t_e(mesh.eLoc(t_e)<0)),:);
    pLocK = mesh.pLoc(mesh.t(tID,:));
    vert0 = mesh.p(mesh.t(tID,:),:);     
    if size(intpt0,1) == 3 % Type 1 Interface Element (3 intersection pts)
        if pLocK(1) == -1 && pLocK(2) == 1 && pLocK(3) == 1 && pLocK(4) == 1 
            locID = [1,2,3,4]; vert = vert0(locID,:); coef = [btM,btP,ahM,ahP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == -1 && pLocK(3) == 1 && pLocK(4) == 1 
            locID = [2,3,4,1]; vert = vert0(locID,:); coef = [btM,btP,ahM,ahP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == 1 && pLocK(3) == -1 && pLocK(4) == 1 
            locID = [3,4,1,2]; vert = vert0(locID,:); coef = [btM,btP,ahM,ahP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == 1 && pLocK(3) == 1 && pLocK(4) == -1 
            locID = [4,1,2,3]; vert = vert0(locID,:); coef = [btM,btP,ahM,ahP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == -1 && pLocK(3) == -1 && pLocK(4) == -1 
            locID = [1,2,3,4]; vert = vert0(locID,:); coef = [btP,btM,ahP,ahM]; plusPC(i) = 1;
        elseif pLocK(1) == -1 && pLocK(2) == 1 && pLocK(3) == -1 && pLocK(4) == -1 
            locID = [2,3,4,1]; vert = vert0(locID,:); coef = [btP,btM,ahP,ahM]; plusPC(i) = 1;
        elseif pLocK(1) == -1 && pLocK(2) == -1 && pLocK(3) == 1 && pLocK(4) == -1 
            locID = [3,4,1,2]; vert = vert0(locID,:); coef = [btP,btM,ahP,ahM]; plusPC(i) = 1;
        elseif pLocK(1) == -1 && pLocK(2) == -1 && pLocK(3) == -1 && pLocK(4) == 1 
            locID = [4,1,2,3]; vert = vert0(locID,:); coef = [btP,btM,ahP,ahM]; plusPC(i) = 1;
        end
        p1 = [vert(1,:);intpt0]; t1 = delaunay(p1);
        p2 = [vert(2:4,:);intpt0];t2 = delaunay(p2);          
    elseif size(intpt0,1) == 4 % Type 2 Interface Element (4 intersection pts)
        if pLocK(1) == -1 && pLocK(2) == -1 && pLocK(3) == 1 && pLocK(4) == 1 
            locID = [1,2,3,4]; vert = vert0(locID,:); coef = [btM,btP,ahM,ahP]; plusPC(i) = 2;
        elseif pLocK(1) == -1 && pLocK(2) == 1 && pLocK(3) == -1 && pLocK(4) == 1 
            locID = [1,3,2,4]; vert = vert0(locID,:); coef = [btM,btP,ahM,ahP]; plusPC(i) = 2;
        elseif pLocK(1) == -1 && pLocK(2) == 1 && pLocK(3) == 1 && pLocK(4) == -1 
            locID = [1,4,2,3]; vert = vert0(locID,:); coef = [btM,btP,ahM,ahP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == -1 && pLocK(3) == -1 && pLocK(4) == 1 
            locID = [1,4,2,3]; vert = vert0(locID,:); coef = [btP,btM,ahP,ahM]; plusPC(i) = 1;
        elseif pLocK(1) == 1 && pLocK(2) == -1 && pLocK(3) == 1 && pLocK(4) == -1 
            locID = [1,3,2,4]; vert = vert0(locID,:); coef = [btP,btM,ahP,ahM]; plusPC(i) = 1;
        elseif pLocK(1) == 1 && pLocK(2) == 1 && pLocK(3) == -1 && pLocK(4) == -1 
            locID = [1,2,3,4]; vert = vert0(locID,:); coef = [btP,btM,ahP,ahM]; plusPC(i) = 1;
        end
        p1 = [vert(1:2,:);intpt0]; t1 = delaunay(p1);
        p2 = [vert(3:4,:);intpt0]; t2 = delaunay(p2);
    end
    
    if tID == 55
        stp1=1;
    end
    
    elocID = 1:6; % reorder the index for edges
    eind = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
    esign = ones(6,1);
    for j = 1:6
        eindnew = locID(eind(j,:));
        if eindnew(2)<eindnew(1)
            esign(j) = -1;
        end
        if min(eindnew) == 1
            elocID(j) = max(eindnew)-1;
        elseif min(eindnew) > 1
            elocID(j) = sum(eindnew)-1;
        end
    end
    if size(intpt0,1) == 3
        intpt = mesh.eIntP(-mesh.eLoc(t_e(elocID(1:3))),:);
    elseif size(intpt0,1) == 4
        intpt = mesh.eIntP(-mesh.eLoc(t_e(elocID(2:5))),:);
    end
    
    
    [BAS1,BAS2] = basIFE3DNed1coef(vert,intpt,coef);     
    X1 = p1(t1(:,1),:); X2 = p1(t1(:,2),:); X3 = p1(t1(:,3),:); X4 = p1(t1(:,4),:);
    A1 = tetraArea(X1,X2,X3,X4);
    badtet1 = find(A1<10^(-14));
    t1(badtet1,:) = [];
    A1(badtet1,:) = [];
    X1 = p1(t1(:,1),:); X2 = p1(t1(:,2),:); X3 = p1(t1(:,3),:); X4 = p1(t1(:,4),:);
    [gx1, gy1, gz1] = gaussPtetra(X1,X2,X3,X4,fem.ng);
    
    X1 = p2(t2(:,1),:); X2 = p2(t2(:,2),:); X3 = p2(t2(:,3),:); X4 = p2(t2(:,4),:);
    A2 = tetraArea(X1,X2,X3,X4);   
    badtet2 = find(A2<10^(-14));
    t2(badtet2,:) = [];
    A2(badtet2,:) = [];
    X1 = p2(t2(:,1),:); X2 = p2(t2(:,2),:); X3 = p2(t2(:,3),:); X4 = p2(t2(:,4),:);
    [gx2, gy2, gz2] = gaussPtetra(X1,X2,X3,X4,fem.ng);
    
    nt1 = size(t1,1); nt2 = size(t2,1);
    bas(id+1:id+nt1,:,elocID(1)) = repmat(BAS1(1,:),nt1,1)*esign(1);
    bas(id+1:id+nt1,:,elocID(2)) = repmat(BAS1(2,:),nt1,1)*esign(2);
    bas(id+1:id+nt1,:,elocID(3)) = repmat(BAS1(3,:),nt1,1)*esign(3);
    bas(id+1:id+nt1,:,elocID(4)) = repmat(BAS1(4,:),nt1,1)*esign(4);
    bas(id+1:id+nt1,:,elocID(5)) = repmat(BAS1(5,:),nt1,1)*esign(5);
    bas(id+1:id+nt1,:,elocID(6)) = repmat(BAS1(6,:),nt1,1)*esign(6);
    bas(id+nt1+1:id+nt1+nt2,:,elocID(1)) = repmat(BAS2(1,:),nt2,1)*esign(1);
    bas(id+nt1+1:id+nt1+nt2,:,elocID(2)) = repmat(BAS2(2,:),nt2,1)*esign(2);
    bas(id+nt1+1:id+nt1+nt2,:,elocID(3)) = repmat(BAS2(3,:),nt2,1)*esign(3);
    bas(id+nt1+1:id+nt1+nt2,:,elocID(4)) = repmat(BAS2(4,:),nt2,1)*esign(4);
    bas(id+nt1+1:id+nt1+nt2,:,elocID(5)) = repmat(BAS2(5,:),nt2,1)*esign(5);
    bas(id+nt1+1:id+nt1+nt2,:,elocID(6)) = repmat(BAS2(6,:),nt2,1)*esign(6);
    gx(id+1:id+nt1+nt2,:) = [gx1;gx2];
    gy(id+1:id+nt1+nt2,:) = [gy1;gy2];
    gz(id+1:id+nt1+nt2,:) = [gz1;gz2];
    A(id+1:id+nt1+nt2,:) = [A1;A2];
    tIntfID(id+1:id+nt1+nt2,:) = tID;
    pLocKnew = pLocK(locID);
    tetID(i,1) = pLocKnew(1);
    tetID(i,2) = pLocKnew(4);
%     if pLocKnew(1) == pLocKnew(4)
%         stp=1;
%     end
%     if i == 85
%         stp=1;
%     end
    
    temp = fem.g2ldof(tID,:);
    temp2 = temp;%(elocID);
    g2ldof(id+1:id+nt1+nt2,:) = repmat(temp2,nt1+nt2,1);
    locIDvec(i,:) = locID;
    elocIDvec(i,:) = elocID;
    quadT(i) = nt1+nt2;
    
    bas1(i,:,elocID(1)) = BAS1(1,:)*esign(1); bas2(i,:,elocID(1)) = BAS2(1,:)*esign(1);
    bas1(i,:,elocID(2)) = BAS1(2,:)*esign(2); bas2(i,:,elocID(2)) = BAS2(2,:)*esign(2);
    bas1(i,:,elocID(3)) = BAS1(3,:)*esign(3); bas2(i,:,elocID(3)) = BAS2(3,:)*esign(3);
    bas1(i,:,elocID(4)) = BAS1(4,:)*esign(4); bas2(i,:,elocID(4)) = BAS2(4,:)*esign(4);
    bas1(i,:,elocID(5)) = BAS1(5,:)*esign(5); bas2(i,:,elocID(5)) = BAS2(5,:)*esign(5);
    bas1(i,:,elocID(6)) = BAS1(6,:)*esign(6); bas2(i,:,elocID(6)) = BAS2(6,:)*esign(6);
    
    id = id+nt1+nt2;
end

eID = find(mesh.eLoc<0);
eIntP = mesh.eIntP;
nge = size(fem.gex,2);
ploc = mesh.pLoc(mesh.e(eID,1)); ploc1 = find(ploc>0);
[gew1tmp,gex1tmp,gey1tmp,gez1tmp] = gaussPedge(mesh.p(mesh.e(eID,1),:),eIntP,nge);
[gew2tmp,gex2tmp,gey2tmp,gez2tmp] = gaussPedge(eIntP,mesh.p(mesh.e(eID,2),:),nge);
gew1 = gew1tmp; gex1 = gex1tmp; gey1 = gey1tmp; gez1 = gez1tmp;
gew2 = gew2tmp; gex2 = gex2tmp; gey2 = gey2tmp; gez2 = gez2tmp;
gew1(ploc1,:) = gew2tmp(ploc1,:); gex1(ploc1,:) = gex2tmp(ploc1,:);
gey1(ploc1,:) = gey2tmp(ploc1,:); gez1(ploc1,:) = gez2tmp(ploc1,:);
gew2(ploc1,:) = gew1tmp(ploc1,:); gex2(ploc1,:) = gex1tmp(ploc1,:);
gey2(ploc1,:) = gey1tmp(ploc1,:); gez2(ploc1,:) = gez1tmp(ploc1,:);

g2ldof(id+1:end,:) = []; bas(id+1:end,:,:) = []; A(id+1:end) = []; tIntfID(id+1:end) = [];
gx(id+1:end,:) = []; gy(id+1:end,:) = []; gz(id+1:end,:) = [];
femI = struct('g2ldof',g2ldof,'bas',bas,'gx',gx,'gy',gy,'gz',gz,'area',A,'gw',fem.gw,...
    'quadT',quadT,'locID',locIDvec,'bas1',bas1,'bas2',bas2,'plusPC',plusPC,...
    'elocIDvec',elocIDvec,'tIntfID',tIntfID,...
    'gew1',gew1,'gex1',gex1,'gey1',gey1,'gez1',gez1,...
    'gew2',gew2,'gex2',gex2,'gey2',gey2,'gez2',gez2,'tetID',tetID);
