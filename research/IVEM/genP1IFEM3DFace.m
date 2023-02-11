function femIF = genP1IFEM3DFace(mesh,fem,femI)
%% Usage: Generate Quadrature Information on Interface Faces Used in PPIFEM
%         Each interface face (triangle) is cut into 3 small triangles
%            A1                A1
%           D  E      or      E  D         equivalent
%          A2  A3            A3   A2
%
%          p = [A1;D;E;A2;A3];
%          t = [1 2 3; 2 3 4; 3 4 5]
%
% femIF.tL --- global index of left element
% femIF.tR --- global index of right element
% femIF.basL --- 3*nfi-4-4 matrix basis functions of left elem on small tri
% femIF.basR --- 3*nfi-4-4 matrix basis functions of right elem on small tri
% femIF.gx --- Gaussian x nodes on small tri
% femIF.gy --- Gaussian y nodes on small tri
% femIF.gz --- Gaussian z nodes on small tri
% femIF.gw --- Gaussian weights on small tri (3 Gaussian pt for linear IFEM)
% femIF.area --- Areas of all small tri
% femIF.normal --- unit normal vector of each small tri

% Last Modified by Xu Zhang 08/07/2020
%%
nfI = -min(mesh.fLoc); % number of interface faces
intFID = find(mesh.fLoc<0);
nfIB = size(find(mesh.f_t(intFID,2)==0),1); % # of interface faces on Boundary
nfII = nfI - nfIB; % # of internal interface faces: two tetra share it L&R
tL = zeros(3*nfII,4); tR = zeros(3*nfII,4);
basL = zeros(3*nfII,4,4); basR = zeros(3*nfII,4,4);
gx = zeros(3*nfII,3); gy = zeros(3*nfII,3); gz = zeros(3*nfII,3);
A = zeros(3*nfII,1); normal = zeros(3*nfII,3);

tB = zeros(3*nfIB,4); basB = zeros(3*nfIB,4,4); 
gxB = zeros(3*nfIB,3); gyB = zeros(3*nfIB,3); gzB = zeros(3*nfIB,3);
AB = zeros(3*nfIB,1); normalB = zeros(3*nfIB,3);

id = 0; idB = 0;
for i = 1:nfI
    %% Form small triangular partition on the interface face
    fID = intFID(i); % face index
    f_e = mesh.f_e(fID,:); % three surrounding edge index
    idE = find(mesh.eLoc(f_e)<0); % find index of two interface edges
    tmp = [mesh.e(f_e(idE(1)),:), mesh.e(f_e(idE(2)),:)];
    nd1 = sum(tmp) - sum(unique(tmp)); % the node of two interface edges
    nd2 = sum(mesh.e(f_e(idE(1)),:)) - nd1;
    nd3 = sum(mesh.e(f_e(idE(2)),:)) - nd1;
    p = [mesh.p(nd1,:); mesh.eIntP(-mesh.eLoc(f_e(idE(1))),:); ...
        mesh.eIntP(-mesh.eLoc(f_e(idE(2))),:); mesh.p(nd2,:); mesh.p(nd3,:)];
    t = [1 2 3; 2 3 4; 3 4 5];
    
    %% gx gy gz on a triangle with 3 internal point, accurate upto pd = 2
    X1 = p(t(:,1),:);  X2 = p(t(:,2),:);   X3 = p(t(:,3),:);
    G = zeros(3,9);  w1 = 2/3;  w2 = 1/6; % see gaussPtri.m
    G(:,[1,4,7]) = w1*X1 + w2*(X2+X3);
    G(:,[2,5,8]) = w1*X2 + w2*(X1+X3);
    G(:,[3,6,9]) = w1*X3 + w2*(X1+X2);
    
    %% triangle area on three-dimension.
%     x1 = X1(:,1);  y1 = X1(:,2);   z1 = X1(:,3);
%     x2 = X2(:,1);  y2 = X2(:,2);   z2 = X2(:,3);
%     x3 = X3(:,1);  y3 = X3(:,2);   z3 = X3(:,3);
%     AT = 1/2*(((x1-x3).*(y2-y1) - (x1-x2).*(y3-y1)).^2 + ...
%         ((y1-y3).*(z2-z1) - (y1-y2).*(z3-z1)).^2 + ...
%         ((z1-z3).*(x2-x1) - (z1-z2).*(x3-x1)).^2).^(1/2);
    AT = TriArea3D(X1,X2,X3);
    
    %% Left and Right Element
    tIDL = mesh.f_t(fID,1); % element index of left element
    tIDR = mesh.f_t(fID,2); % element index of right element
    if tIDR > 0 % Internal Face 
        tIDLi = -mesh.tLoc(tIDL); % intf elem index of left element
        tIDRi = -mesh.tLoc(tIDR); % intf elem index of right element
        
        %% Determine piece
        nd1ID = mesh.pLoc(nd1); tLp = femI.plusPC(tIDLi); tRp = femI.plusPC(tIDRi);
        if (nd1ID < 0 && tLp == 1) || (nd1ID > 0 && tLp == 2)
            basL(id+1,:,:) = femI.bas2(tIDLi,:,:);
            basL(id+2,:,:) = femI.bas1(tIDLi,:,:);
            basL(id+3,:,:) = femI.bas1(tIDLi,:,:);
        elseif (nd1ID < 0 && tLp == 2) || (nd1ID > 0 && tLp == 1)
            basL(id+1,:,:) = femI.bas1(tIDLi,:,:);
            basL(id+2,:,:) = femI.bas2(tIDLi,:,:);
            basL(id+3,:,:) = femI.bas2(tIDLi,:,:);
        end
        if (nd1ID < 0 && tRp == 1) || (nd1ID > 0 && tRp == 2)
            basR(id+1,:,:) = femI.bas2(tIDRi,:,:);
            basR(id+2,:,:) = femI.bas1(tIDRi,:,:);
            basR(id+3,:,:) = femI.bas1(tIDRi,:,:);
        elseif (nd1ID < 0 && tRp == 2) || (nd1ID > 0 && tRp == 1)
            basR(id+1,:,:) = femI.bas1(tIDRi,:,:);
            basR(id+2,:,:) = femI.bas2(tIDRi,:,:);
            basR(id+3,:,:) = femI.bas2(tIDRi,:,:);
        end
        
        gx(id+1:id+3,:) = G(:,1:3);
        gy(id+1:id+3,:) = G(:,4:6);
        gz(id+1:id+3,:) = G(:,7:9);
        
        A(id+1:id+3,:) = AT;
        normal(id+1:id+3,:) = repmat(mesh.f_norm(fID,:),3,1);
        
        %% tL and tR, use locID, b/c index on interface cell is different
        temp = fem.t(tIDL,:);
        temp1 = temp(femI.locID(tIDLi,:));
        tL(id+1:id+3,:) = repmat(temp1,3,1);
        
        temp = fem.t(tIDR,:);
        temp2 = temp(femI.locID(tIDRi,:));
        tR(id+1:id+3,:) = repmat(temp2,3,1);
        id = id+3;
        
    elseif tIDR == 0 % Boundary Face
        tIDLi = -mesh.tLoc(tIDL); % intf elem index of left element
        
        %% Determine piece: only one element.
        nd1ID = mesh.pLoc(nd1); tLp = femI.plusPC(tIDLi); 
        if (nd1ID < 0 && tLp == 1) || (nd1ID > 0 && tLp == 2)
            basB(idB+1,:,:) = femI.bas2(tIDLi,:,:);
            basB(idB+2,:,:) = femI.bas1(tIDLi,:,:);
            basB(idB+3,:,:) = femI.bas1(tIDLi,:,:);
        elseif (nd1ID < 0 && tLp == 2) || (nd1ID > 0 && tLp == 1)
            basB(idB+1,:,:) = femI.bas1(tIDLi,:,:);
            basB(idB+2,:,:) = femI.bas2(tIDLi,:,:);
            basB(idB+3,:,:) = femI.bas2(tIDLi,:,:);
        end
        
        gxB(idB+1:idB+3,:) = G(:,1:3);
        gyB(idB+1:idB+3,:) = G(:,4:6);
        gzB(idB+1:idB+3,:) = G(:,7:9);
        
        AB(idB+1:idB+3,:) = AT;
        normalB(idB+1:idB+3,:) = repmat(mesh.f_norm(fID,:),3,1);
        
        %% tB use locID, b/c index on interface cell is different
        temp = fem.t(tIDL,:);
        temp1 = temp(femI.locID(tIDLi,:));
        tB(idB+1:idB+3,:) = repmat(temp1,3,1);
        idB = idB+3;
    end
end

femIF = struct('tL',tL,'tR',tR,'tB',tB,'basL',basL,'basR',basR,'basB',basB, ...
    'gx',gx,'gy',gy,'gz',gz,'gxB',gxB,'gyB',gyB,'gzB',gzB,'area',A,'areaB',AB,...
    'gw',[1/3;1/3;1/3],'normal',normal,'normalB',normalB);
