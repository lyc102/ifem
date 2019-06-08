function surfacedata = torussurface
R = 4;
r = 1;
Nu = 20;
Nv = 10;
surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @gradient, 'project', @project, 'Hessian', @H, 'tanproperator',...
    @tanproperator, 'initmesh', @initmesh, 'meancurvature', @meancurvature);


%---------------------------------------------------

function [node,elem]= initmesh

M = load('meshdata/torus4770.mat');
node = M.node;
elem = M.elem;



% u=(0:2*pi/Nu:2*pi)';
% v=(0:2*pi/Nv:2*pi)';
% elem=zeros(2*Nu*Nv,3);
% count =0;
% for i=1:Nu
%     for j=1:Nv
%         if (i<Nu && j<Nv)
%             count = count + 1;
%             elem(count,:)=[(i-1)*Nv+j+1,(i-1)*Nv+j,i*Nv+j+1];
%             count = count + 1;
%             elem(count,:)=[i*Nv+j, i*Nv+j+1, (i-1)*Nv+j];
%         elseif (i==Nu &&j<Nv)
%             count = count + 1;
%             elem(count,:)=[(i-1)*Nv+j+1,(i-1)*Nv+j, j+1];
%             count = count + 1;
%             elem(count,:)=[j,j+1,(i-1)*Nv+j];
%         elseif (i<Nu && j==Nv)
%             count = count + 1;            
%             elem(count,:)=[(i-1)*Nv+1,(i-1)*Nv+j, i*Nv+1];
%             count = count + 1;
%             elem(count,:)=[i*Nv+j, i*Nv+1 ,(i-1)*Nv+j];
%         else
%             count = count + 1;
%             elem(count,:)=[(i-1)*Nv+1,(i-1)*Nv+j,1];
%             count = count + 1;
%             elem(count,:)=[j,1,(i-1)*Nv+j];
%         end
%     end
% end
% 
% node=zeros(Nu*Nv,3);
% count = 0;
% for i=1:Nu
%     for j=1:Nv
%         count = count + 1;
%         node(count,:)=[(R+r*cos(v(j)))*cos(u(i)),(R+r*cos(v(j)))*sin(u(i)),r*sin(v(j))];
%     end
% end


function [z, dist] = project(p)
% projection function 
d = phi(p);
z = p - d*ones(1,3).*gradient(p);
if nargout == 2
    dist = d;
end




function z = phi(p)

z = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2+16-8*sqrt(p(:,1).^2+p(:,2).^2))-1;


function n = gradient(p)

s1=sqrt(p(:,1).^2+p(:,2).^2);
s2=sqrt(s1.^2 + p(:,3).^2+16-8*s1);
n1 = -(4.0 - s1).* p(:,1)./(s1.*s2);
n2 = -(4.0 - s1).* p(:,2)./(s1.*s2);
n3 = p(:,3)./s2;
n = [n1,n2,n3];


function H = Hessian(p)

H=zeros(3,3,size(p,1));

s = sqrt(p(:,1).^2+p(:,2).^2);

d11 = -(p(:,1) - 4 .* p(:,1) ./ s) .^ 2 + 1 + 4 .* p(:,1) .^ 2 ./ s .^ 3 - 4 ./ s;

d12 = p(:,1) .* p(:,2) .* (1 - 4 ./ s) .^ 2 + 4 .* p(:,1) .* p(:,2) ./ s .^ 3;

d13= -(p(:,1) - 4 .* p(:,1) ./ s) .* p(:,3);

d22 = (p(:,2) - 4 .* p(:,2) ./ s) .^ 2 + 1 + 4 .* p(:,2) .^ 2 ./ s .^ 3 - 4 ./ s;

d23 = (p(:,2) - 4 .* p(:,2) ./ s) .* p(:,3);

d33 = -p(:,3) .^ 2 + 1;

for i=1:size(p,1)
    H(:,:,i)=[d11(i),d12(i),d13(i);d12(i),d22(i),d23(i);d13(i),d23(i),d33(i)];
end

function mc = meancurvature(p)
% H = Hessian(p);
% N = length(p);
% mc = zeros(N,1);
% for i = 1:N
%     mc(i) = trace(H(:,:,i))/2;
% end

x = p(:,1);
y = p(:,2);
z = p(:,3);
s = sqrt(x.^2+y.^2);
s1 = x.^2+y.^2+z.^2+16;
s2 = (s1+16).*s-2*z.^2 - 10*x.^2 - 10*y.^2 - 32;
s3 =s.*(s1-8*s).^(3/2);
mc = s2./s3;

function z = tanproperator(p)
z=zeros(3,3,size(p,1));
p = project(p);
n = gradient(p);
for i=1:size(p,1)
   z(:,:,i)=n(i,:)'*n(i,:); 
end

