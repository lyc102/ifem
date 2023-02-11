elemID = 1;
[Dlambda,volume,elemSign] = gradbasis3(node,mesh.t(elemID,:));
elem2face = [1,2,3,4]; 
localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
elem2faceSignLoc = [1,1,1,1];
MatLoc = sparse(4,4);
for i = 1:4
    for j = i:4
        % local to global index map
        ii = double(elem2face(:,i));
        jj = double(elem2face(:,j));
        i1 = localFace(i,1); i2 = localFace(i,2); i3 = localFace(i,3);
        j1 = localFace(j,1); j2 = localFace(j,2); j3 = localFace(j,3);
        % computation of mass matrix --- (phi_i, phi_j)
        Mij = 1/20*volume*4.*( ...
            (1+(i1==j1))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
            mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
            +(1+(i1==j2))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
            mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
            +(1+(i1==j3))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
            mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)...
            +(1+(i2==j1))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
            mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
            +(1+(i2==j2))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
            mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
            +(1+(i2==j3))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
            mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)...
            +(1+(i3==j1))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
            mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
            +(1+(i3==j2))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
            mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
            +(1+(i3==j3))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
            mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)).*...
            elem2faceSignLoc(:,i).*elem2faceSignLoc(:,j);
            %elem2faceArea(:,i).^(-1).*elem2faceArea(:,j).^(-1).*...
        if (j==i)
            MatLoc = MatLoc + sparse(ii,jj,Mij,4,4);
        else
            MatLoc = MatLoc + sparse([ii;jj],[jj;ii],[Mij; Mij],4,4);
        end
    end
end

MatLoc = full(MatLoc);


MatInd = 1;
ntID = elemID;
AN = fem.area(ntID); 
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); gw = fem.gw;
XN = zeros(36, 1);

coefN = feval(pde.A,gxN,gyN,gzN);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);

for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 1);
    Ibasy{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 2);
    Ibasz{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 3);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 1);
    Jbasy{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 2);
    Jbasz{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 3);
end

IN = reshape(repmat(1:6,6,1),36,1);
JN = repmat(reshape(1:6,6,1),6,1);
ind = 0;
for i = 1:dof1
    for j = 1:dof2
        XN(ind+1:ind+1) = AN.*(sum(((Ibasx{i}.*(coefN.*Jbasx{j})).*gw'),2) + ...
            sum(((Ibasy{i}.*(coefN.*Jbasy{j})).*gw'),2) + ...
            sum(((Ibasz{i}.*(coefN.*Jbasz{j})).*gw'),2));
        ind = ind + 1;
    end
end
SN = sparse(IN,JN,XN,6,6);
SN = full(SN);

ICLoc = [-1, 1, 0,-1, 0, 0;...
          1, 0,-1, 0, 1, 0;...
          0,-1, 1, 0, 0,-1;...
          0, 0, 0, 1,-1,1];

% faceLoc = [1,2,3;1,2,4;1,3,4;2,3,4];
% checkid = zeros(size(mesh.t,1),1);
% for i = 1:size(mesh.t,1)
%     for j = 1:4
%         v1=mesh.f(mesh.t_f(i,j),:);
%         v2=mesh.t(i,faceLoc(j,:));
%         if sort(v1)~=sort(v2)
%             checkid(i)=1;
%         end
%     end
% end


