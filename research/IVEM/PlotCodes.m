error1 = Error(1,:);
error2 = Error(2,:);
Intpt =10.^(-1:-1:-6);

semilogx(Intpt,error1,'r-o')
hold on
semilogx(Intpt,error2,'r-*')

lgd = legend('L^2 error','H(curl) error');
lgd.FontSize = 14;


xlim([10^(-7) 1]) 
ylim([0.5 1.2])

semilogx(Intpt,ConNum,'r-s')
ylim([7*10^4 10*10^4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lev = 8;
N = [10:10:10*lev];
h = 2./[10:10:10*lev]';
load('PGIFE_error_am1ap1000_bm1bp1000.mat')
% load('PGIFE_error_am1ap1000_bm1bp1000.mat')
% error2 = error;
r10 = regress(log(error(:,2)),[ones(length(h),1),log(h)])
r11 = regress(log(error(:,3)),[ones(length(h),1),log(h)])
% r20 = regress(log(error2(:,2)),[ones(length(h),1),log(h)])
% r21 = regress(log(error2(:,3)),[ones(length(h),1),log(h)])
expec_err0=zeros(1,lev);
expec_err0(lev)=error(end,2)/2;
for j=1:lev-1
    expec_err0(j)=(N(lev)/N(j))*expec_err0(lev);
end
loglog(N,error(:,2),'r-o')
hold on
loglog(N,error(:,3),'b-*')
loglog(N, expec_err0,'k--');
lgd = legend('L^2 error','H(curl) error');
lgd.FontSize = 14;
text(70,0.6,'\fontname{Courier Oblique} h^{ 1.09}','FontSize',20,'color','r')
text(70,2.5,'\fontname{Courier Oblique} h^{ 1.14}','FontSize',20,'color','b')

axis([8 100 10^(-1) 10^(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lev = 16;
N = [10:10:10*lev];
h = 2./[10:10:10*lev]';
load('H1_error_bm1bp100_rpi4.mat')
% load('PGIFE_error_am1ap1000_bm1bp1000.mat')
% error2 = error;
stlev = 6;
r10 = regress(log(error(stlev:lev,2)),[ones(length(h(stlev:lev)),1),log(h(stlev:lev))])
r11 = regress(log(error(stlev:lev,3)),[ones(length(h(stlev:lev)),1),log(h(stlev:lev))])
r12 = regress(log(error(stlev:lev,4)),[ones(length(h(stlev:lev)),1),log(h(stlev:lev))])
% r20 = regress(log(error2(:,2)),[ones(length(h),1),log(h)])
% r21 = regress(log(error2(:,3)),[ones(length(h),1),log(h)])
expec_err0=zeros(1,lev);
expec_err0(lev)=error(lev,2);
expec_err1=zeros(1,lev);
expec_err1(lev)=error(lev,3);
expec_err2=zeros(1,lev);
expec_err2(lev)=error(lev,4);
for j=1:lev-1
    expec_err0(j)=(N(lev)/N(j))^2*expec_err0(lev);
    expec_err1(j)=(N(lev)/N(j))^2*expec_err1(lev);
    expec_err2(j)=(N(lev)/N(j))*expec_err2(lev);
end
loglog(N,error(1:lev,2),'r-o')
hold on
loglog(N,error(1:lev,3),'b-*')
loglog(N,error(1:lev,4),'m-s')
loglog(N, expec_err0,'k--');
loglog(N, expec_err1,'k--');
loglog(N, expec_err2,'k--');
lgd = legend('L^{\infty} error','L^2 error','H^1 error');
lgd.FontSize = 14;
text(170,2.5*10^(-4),'\fontname{Courier Oblique} h^{ 1.90}','FontSize',20,'color','r')
text(170,1.2*10^(-4),'\fontname{Courier Oblique} h^{ 1.99}','FontSize',20,'color','b')
text(170,0.019,'\fontname{Courier Oblique} h^{ 0.97}','FontSize',20,'color','m')

axis([8 300 5*10^(-5) 0.4])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot infsup constants
%h = 2./[2,4,8,16,24,32,40,48];
h = 10:10:70;
plot(h,sqrt(EigAll1),'r-o','MarkerSize',10,'MarkerEdgeColor','red')
hold on
plot(h,sqrt(EigAll2),'b-*','MarkerSize',10,'MarkerEdgeColor','blue')
plot(h,sqrt(EigAll3),'-^','MarkerSize',10,'Color',[0.5 0 0.8])

axis([8,70,0,0.3])

le = legend('\beta^+=200,\alpha^+=100','\beta^+=200,\alpha^+=1000','\beta^+=2000,\alpha^+=100');
le.FontSize = 15;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot interface surface with triangles
mpt1 = (node(tface(:,1),:) + node(tface(:,2),:) + node(tface(:,3),:))/3;
ptId1 = mpt1(:,1)<=0;
trimesh(tface(ptId1,:),node(:,1),node(:,2),node(:,3),'FaceAlpha',0,'EdgeColor','b')
hold on
mpt2 = (node(iface(:,1),:) + node(iface(:,2),:) + node(iface(:,3),:))/3;
ptId2 = mpt2(:,1)>=0;
trimesh(iface(ptId2,:),node(:,1),node(:,2),node(:,3),'FaceAlpha',0,'EdgeColor','r')
domain = [-1,1,-1,1,-1,1];
axis(domain)
axis equal



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = @(x,y,z)(((x+0.3).^2+y.^2).^(1/2)-pi/5).^2+z.^2-0.2^2;
f2 = @(x,y,z)(((x-0.3).^2+z.^2).^(1/2)-pi/5).^2+y.^2-0.2^2;
f = @(x,y,z) f1(x,y,z).*f2(x,y,z);
fimplicit3(f)
axis equal
domain = [-1.3,1.3,-1.3,1.3,-1.3,1.3];
axis(domain)







