function surfacedata = spheresurface

surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @gradient, 'project', @project, 'Hessian', @Hessian, ...
    'tanproperator', @tanproperator, 'initmesh',@initmesh);


function [node,elem] = initmesh
%Produce a icosahedron mesh.
t=(sqrt(5)-1)/2;
node=[0, 1, t
      0, 1,-t
      1, t, 0
      1,-t, 0
      0,-1,-t
      0,-1, t
      t, 0, 1
      -t,0, 1
      t, 0,-1
      -t,0,-1
      -1,t, 0
      -1,-t,0];
 elem=[ 6, 2,0
        3, 2, 6
        5, 3, 6
        5, 6, 7
        6, 0, 7
        3, 8, 2
        2, 8, 1
        2, 1, 0
        0, 1, 10
        1, 9, 10
        8, 9, 1
        4, 8, 3
        4, 3, 5
        4, 5, 11
        7, 10, 11
        0, 10, 7
        4, 11, 9
        8, 4, 9
        5, 7, 11
        10, 9, 11]+1;
 node= project(node);
 
 

function [z,dist] = project(p)
% projection function 
d = phi(p);
z = p - d*ones(1,3).*gradient(p);
if nargout == 2
   dist = d; 
end



function z = phi(p)
% level set function \phi
z = p(:,1).^2+p(:,2).^2+p(:,3).^2-0.75^2;


function n = gradient(p)
% the gradient function of \phi
L = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
n = p ./(L*ones(1,3));



function H = Hessian(p)
% Hessian matrix
H = zeros(3,3,size(p,1));
L = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
L3=L.^3;
d11 = 1./L-p(:,1).^2./L3;
d12 = -p(:,1).*p(:,2)./L3;
d13 = -p(:,1).*p(:,3)./L3;
d22 = 1./L-p(:,2).^2./L3;
d23 = -p(:,2).*p(:,3)./L3;
d33 = 1./L-p(:,3).^2./L3;
for i=1:size(p,1)
    H(:,:,i)=[d11(i),d12(i),d13(i);d12(i),d22(i),d23(i);d13(i),d23(i),d33(i)];
end

function z = tanproperator(p)

z=zeros(3,3,size(p,1));
n = gradient(p);
for i=1:size(p,1)
   z(:,:,i)=n(i,:)'*n(i,:); 
end