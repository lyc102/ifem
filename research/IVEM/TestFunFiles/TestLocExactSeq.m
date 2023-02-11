p = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
t = [1,2,3,4];
bas = bas3DP1(p,t);
bas = bas3D_ned1(bas);
S = zeros(6,6);
feEvalBas = @EvalNed1Bas3D;
gw = gaussWtetra(4);
[gx,gy,gz] = gaussPtetra(p(1,:),p(2,:),p(3,:),p(4,:),4);

for i = 1:6
    for j = 1:6
        
        bix = feEvalBas(bas, ':', gx, gy, gz, i, 1, 1);
        biy = feEvalBas(bas, ':', gx, gy, gz, i, 1, 2);
        biz = feEvalBas(bas, ':', gx, gy, gz, i, 1, 3);
        
        bjx = feEvalBas(bas, ':', gx, gy, gz, j, 1, 1);
        bjy = feEvalBas(bas, ':', gx, gy, gz, j, 1, 2);
        bjz = feEvalBas(bas, ':', gx, gy, gz, j, 1, 3);
        
        S(i,j) = sum((bix.*bjx + biy.*bjy + biz.*bjz).*gw);
        
    end
end

G = [1,-1,0,0; 1,0,-1,0; 1,0,0,-1; 0,1,-1,0; 0,1,0,-1; 0,0,1,-1];
%l = sqrt([1,1,1,sqrt(2),sqrt(2),sqrt(2)]);
G = G;%.*(l'*ones(1,4));
E = G'*S*G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test the local exact sequence on interface elements
% # cutting points = 3
vert = p;
intpt = [1/2,0,0; 0,1/2,0; 0,0,1/2];
coef = [1,1,1,1];
p1 = [p(1,:);intpt];
p2 = [p(2,:);intpt];
t1 = delaunay(p1);
t2 = delaunay(p2);
X1 = p1(t1(:,1),:); X2 = p1(t1(:,2),:); X3 = p1(t1(:,3),:); X4 = p1(t1(:,4),:);
[gxloc1, gyloc1, gzloc1] = gaussPtetra(X1,X2,X3,X4,fem.ng);
Aloc1 = tetraArea(X1,X2,X3,X4);
X1 = p2(t2(:,1),:); X2 = p2(t2(:,2),:); X3 = p2(t2(:,3),:); X4 = p2(t2(:,4),:);
[gxloc2, gyloc2, gzloc2] = gaussPtetra(X1,X2,X3,X4,fem.ng);
Aloc2 = tetraArea(X1,X2,X3,X4);
OrietCorrect = 1;
sign_correct = ones(12,1);
[BAScurl1,BAScurl2,KCurlstiff,SCurl,BASu1,BASu2,KUstiff,SU] =...
            basVIFE3DNed1coef1(vert,intpt,coef,Aloc1,Aloc2,p1,p2,t1,t2,...
            OrietCorrect,sign_correct);
        
EdgeEndID = [1,1,1,2,2,3,5,6,7,5,5,6;
    5,6,7,3,4,4,2,3,4,6,7,7]';
G = zeros(size(EdgeEndID,1));
for j = 1:size(EdgeEndID,1)
    G(j,EdgeEndID(j,:)) = [1,-1];
end
E = G'*KCurlstiff*G;
% v = ones(1,size(EdgeEndID,1));
% for i = 0:2^(size(EdgeEndID,1))
%     ii = de2bi(i);
%     v(find(ii==0)) = -1;
%     for j = 1:size(EdgeEndID,1)
%         G(j,EdgeEndID(j,:)) = [1,-1]*v(j);
%     end
%     E = G'*KCurlstiff*G;
%     if max(max(abs(E)))<10^(-10)
%         stp = 1
%     end
% end

edge = [1,1,1,2,2,3,5,6,7,5,5,6;
        5,6,7,3,4,4,2,3,4,6,7,7]';
NE = size(edge,1);
N = 7;
i1 = (1:NE)'; j1 = double(edge(:,1)); s1 = ones(NE,1);
i2 = (1:NE)'; j2 = double(edge(:,2)); s2 = -s1;
G = sparse([i1(:);i2(:)],...
    [j1(:);j2(:)],...
    [s1(:);s2(:)],NE,N);
for j = 1:NE
    G(j,:) = G(j,:)*sign_correct(j);
end

             
             
             