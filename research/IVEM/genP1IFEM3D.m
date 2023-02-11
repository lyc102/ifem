function femI = genP1IFEM3D(mesh,fem,btM,btP)
%% Usage: Generate Quadrature Information on Interface Elements
%         Each interface tetrahedron is cut into 4 to 9 small tetrahedra 
%
% femI.t --- global index on small tetrahedra
% femI.bas --- nt-4-4 matrix stores basis functions on all small tetra   
% femI.gx --- Gaussian x nodes
% femI.gy --- Gaussian y nodes
% femI.gz --- Gaussian z nodes
% femI.gw --- Gaussian weights
% femI.area --- Areas of all small tetra
% femI.quadT --- number of small tetras of each interface element
% femI.locID --- order of local index, to match with the basis function

% ***** The following info only used for partial penalty IFEM *****

% femI.bas1 --- nti-4-4 basis function on 1st piece of interface element
% femI.bas2 --- nti-4-4 basis function on 2nd piece of interface element
% femI.plusPC --- 1 or 2 denotes which piece corresponds to btP.

% Last Modified by Xu Zhang 08/02/2020

%% 
ntI = -min(mesh.tLoc);
t = zeros(8*ntI,4); bas = zeros(8*ntI,4,4); A = zeros(8*ntI,1);
gx = zeros(8*ntI,4); gy = zeros(8*ntI,4); gz = zeros(8*ntI,4);
QelemID = zeros(8*ntI,1);
quadT = zeros(ntI,1); locIDvec = zeros(ntI,4);
bas1 = zeros(ntI,4,4); bas2 = zeros(ntI,4,4); plusPC = zeros(ntI,1);
intID = find(mesh.tLoc<0);
id = 0;
for i = 1:ntI
    tID = intID(i);
    t_e = mesh.t_e(tID,:);
    intpt = mesh.eIntP(-mesh.eLoc(t_e(mesh.eLoc(t_e)<0)),:);
    pLocK = mesh.pLoc(mesh.t(tID,:));
    vert0 = mesh.p(mesh.t(tID,:),:);     
    if size(intpt,1) == 3 % Type 1 Interface Element (3 intersection pts)
        if pLocK(1) == -1 && pLocK(2) == 1 && pLocK(3) == 1 && pLocK(4) == 1 
            locID = [1,2,3,4]; vert = vert0(locID,:); coef = [btM,btP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == -1 && pLocK(3) == 1 && pLocK(4) == 1 
            locID = [2,3,4,1]; vert = vert0(locID,:); coef = [btM,btP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == 1 && pLocK(3) == -1 && pLocK(4) == 1 
            locID = [3,4,1,2]; vert = vert0(locID,:); coef = [btM,btP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == 1 && pLocK(3) == 1 && pLocK(4) == -1 
            locID = [4,1,2,3]; vert = vert0(locID,:); coef = [btM,btP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == -1 && pLocK(3) == -1 && pLocK(4) == -1 
            locID = [1,2,3,4]; vert = vert0(locID,:); coef = [btP,btM]; plusPC(i) = 1;
        elseif pLocK(1) == -1 && pLocK(2) == 1 && pLocK(3) == -1 && pLocK(4) == -1 
            locID = [2,3,4,1]; vert = vert0(locID,:); coef = [btP,btM]; plusPC(i) = 1;
        elseif pLocK(1) == -1 && pLocK(2) == -1 && pLocK(3) == 1 && pLocK(4) == -1 
            locID = [3,4,1,2]; vert = vert0(locID,:); coef = [btP,btM]; plusPC(i) = 1;
        elseif pLocK(1) == -1 && pLocK(2) == -1 && pLocK(3) == -1 && pLocK(4) == 1 
            locID = [4,1,2,3]; vert = vert0(locID,:); coef = [btP,btM]; plusPC(i) = 1;
        end
        p1 = [vert(1,:);intpt]; t1 = delaunay(p1);
        p2 = [vert(2:4,:);intpt];t2 = delaunay(p2);        
    elseif size(intpt,1) == 4 % Type 2 Interface Element (4 intersection pts)
        if pLocK(1) == -1 && pLocK(2) == -1 && pLocK(3) == 1 && pLocK(4) == 1 
            locID = [1,2,3,4]; vert = vert0(locID,:); coef = [btM,btP]; plusPC(i) = 2;
        elseif pLocK(1) == -1 && pLocK(2) == 1 && pLocK(3) == -1 && pLocK(4) == 1 
            locID = [1,3,2,4]; vert = vert0(locID,:); coef = [btM,btP]; plusPC(i) = 2;
        elseif pLocK(1) == -1 && pLocK(2) == 1 && pLocK(3) == 1 && pLocK(4) == -1 
            locID = [1,4,2,3]; vert = vert0(locID,:); coef = [btM,btP]; plusPC(i) = 2;
        elseif pLocK(1) == 1 && pLocK(2) == -1 && pLocK(3) == -1 && pLocK(4) == 1 
            locID = [1,4,2,3]; vert = vert0(locID,:); coef = [btP,btM]; plusPC(i) = 1;
        elseif pLocK(1) == 1 && pLocK(2) == -1 && pLocK(3) == 1 && pLocK(4) == -1 
            locID = [1,3,2,4]; vert = vert0(locID,:); coef = [btP,btM]; plusPC(i) = 1;
        elseif pLocK(1) == 1 && pLocK(2) == 1 && pLocK(3) == -1 && pLocK(4) == -1 
            locID = [1,2,3,4]; vert = vert0(locID,:); coef = [btP,btM]; plusPC(i) = 1;
        end
        p1 = [vert(1:2,:);intpt]; t1 = delaunay(p1);
        p2 = [vert(3:4,:);intpt]; t2 = delaunay(p2);
    else
        stp = 1;
    end
    [BAS1,BAS2] = basIFE3DP1coef(vert,intpt,coef);     
    X1 = p1(t1(:,1),:); X2 = p1(t1(:,2),:); X3 = p1(t1(:,3),:); X4 = p1(t1(:,4),:);
    [gx1, gy1, gz1] = gaussPtetra(X1,X2,X3,X4,fem.ng);
    A1 = tetraArea(X1,X2,X3,X4);
    
    X1 = p2(t2(:,1),:); X2 = p2(t2(:,2),:); X3 = p2(t2(:,3),:); X4 = p2(t2(:,4),:);
    [gx2, gy2, gz2] = gaussPtetra(X1,X2,X3,X4,fem.ng);
    A2 = tetraArea(X1,X2,X3,X4);    
    
    nt1 = size(t1,1); nt2 = size(t2,1);
    bas(id+1:id+nt1,:,1) = repmat(BAS1(1,:),nt1,1);
    bas(id+1:id+nt1,:,2) = repmat(BAS1(2,:),nt1,1);
    bas(id+1:id+nt1,:,3) = repmat(BAS1(3,:),nt1,1);
    bas(id+1:id+nt1,:,4) = repmat(BAS1(4,:),nt1,1);
    bas(id+nt1+1:id+nt1+nt2,:,1) = repmat(BAS2(1,:),nt2,1);
    bas(id+nt1+1:id+nt1+nt2,:,2) = repmat(BAS2(2,:),nt2,1);
    bas(id+nt1+1:id+nt1+nt2,:,3) = repmat(BAS2(3,:),nt2,1);
    bas(id+nt1+1:id+nt1+nt2,:,4) = repmat(BAS2(4,:),nt2,1);
    gx(id+1:id+nt1+nt2,:) = [gx1;gx2];
    gy(id+1:id+nt1+nt2,:) = [gy1;gy2];
    gz(id+1:id+nt1+nt2,:) = [gz1;gz2];
    A(id+1:id+nt1+nt2,:) = [A1;A2];
    QelemID(id+1:id+nt1+nt2) = tID;
    
    temp = fem.t(tID,:);
    temp2 = temp(locID);
    t(id+1:id+nt1+nt2,:) = repmat(temp2,nt1+nt2,1);
    locIDvec(i,:) = locID;
    quadT(i) = nt1+nt2;
    
    bas1(i,:,1) = BAS1(1,:); bas2(i,:,1) = BAS2(1,:);
    bas1(i,:,2) = BAS1(2,:); bas2(i,:,2) = BAS2(2,:);
    bas1(i,:,3) = BAS1(3,:); bas2(i,:,3) = BAS2(3,:);
    bas1(i,:,4) = BAS1(4,:); bas2(i,:,4) = BAS2(4,:);
    
    id = id+nt1+nt2;
end
t(id+1:end,:) = []; bas(id+1:end,:,:) = []; A(id+1:end) = [];
gx(id+1:end,:) = []; gy(id+1:end,:) = []; gz(id+1:end,:) = [];
QelemID(id+1:end,:) = [];
femI = struct('t',t,'bas',bas,'gx',gx,'gy',gy,'gz',gz,'area',A,'gw',fem.gw,...
    'quadT',quadT,'locID',locIDvec,'bas1',bas1,'bas2',bas2,'plusPC',plusPC,...
    'QelemID',QelemID);
