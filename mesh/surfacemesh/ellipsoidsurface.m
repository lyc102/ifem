function surfacedata = ellipsoidsurface

a = 9;
b = 3;
c = 1;

a2 = a*a;
a4 = a2*a2;
b2 = b*b;
b4 = b2*b2;
c2 = c*c;
c4 = c2*c2;

surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'initmesh',@initmesh, 'meancurvature', @meancurvature, ...
    'u',@u, 'f',@f,'gra_u',@gra_u);


function [node, elem] = initmesh
M = load('meshdata/ellipsoid3394.mat');
node = M.node;
elem = M.elem;
end


function z = phi(p)
% level set function

x = p(:,1); y = p(:,2); z = p(:,3);
z = x.^2/a2 + y.^2/b2 + z.^2/c2 - 1.0;
end


function n = unitoutnormal(p)
% 
n = gradient(p);
l = sqrt(sum(n.^2,2));
n = n./repmat(l,1,3);
end

function n = gradient(p)
% gradient of phi

x = p(:,1); y = p(:,2); z = p(:,3);
n = [2*x/a2, 2*y/b2, 2*z/c2];
end

function [node,dist] = project(p)
% projection function 

s = sign(phi(p));

node = p;

normalAtNode = gradient(node);
valueAtNode = phi(node);
node = node - valueAtNode*ones(1,3).*normalAtNode./(dot(normalAtNode,normalAtNode,2)*ones(1,3));

vector = (-s)*ones(1,3).*(node - p);

d = s.*sqrt(dot(vector,vector,2));

normalAtNode = gradient(node);

node = p - d*ones(1,3).*normalAtNode./(sqrt(dot(normalAtNode,normalAtNode,2))*ones(1,3));

valueAtNode = phi(node);

normalAtNode = gradient(node);

vector = (-s)*ones(1,3).*(node - p);

d = s.*sqrt(dot(vector,vector,2));


e1 = normalAtNode./(sqrt(dot(normalAtNode, normalAtNode,2))*ones(1,3))-vector./(sqrt(dot(vector,vector,2))*ones(1,3));
error=sqrt(valueAtNode.^2./(dot(normalAtNode, normalAtNode,2))+dot(e1,e1,2));

k=1;
while max(abs(error)) > 1e-6 && k<200
    
    k=k+1;
    
    node = node - valueAtNode*ones(1,3).*normalAtNode./(dot(normalAtNode,normalAtNode,2)*ones(1,3));
    
    vector = -s*ones(1,3).*(node - p);
    d = s.*sqrt(dot(vector,vector,2));
    normalAtNode = gradient(node);
    
    node = p - d*ones(1,3).*normalAtNode./(sqrt(dot(normalAtNode,normalAtNode,2))*ones(1,3));
    
    valueAtNode = phi(node);
    
    normalAtNode = gradient(node);
    
    vector = (-s)*ones(1,3).*(node - p);
    
    d = s.*sqrt(dot(vector,vector,2));
    e1 = normalAtNode./(sqrt(dot(normalAtNode, normalAtNode,2))*ones(1,3))-vector./(sqrt(dot(vector,vector,2))*ones(1,3));
    error=sqrt(valueAtNode.^2./(dot(normalAtNode, normalAtNode,2))+dot(e1,e1,2));   
end

if nargout == 2
    dist = d;
end
end


function mc = meancurvature(p)

x = p(:,1); y = p(:,2); z = p(:,3);
x2 = x.^2; y2 = y.^2; z2 = z.^2;


s1 = c4*a2*y2 + b4*a2*z2 + c4*b2*x2 + a4*b2*z2 + b4*c2*x2 + a4*c2*y2;
s2 = b4*c4*x2 + a4*c4*y2 + a4*b4*z2;
s3 = x2/a4+y2/b4+z2/c4;

mc = 0.5* s1./(s2.*sqrt(s3));
end


function z = f(p)
p = project(p);
x = p(:,1); y = p(:,2); z = p(:,3);
x2 = x.^2; y2 = y.^2; z2 = z.^2;
d1 = c4*y2 + b4*z2;
d2 = c2*y2 + b2*z2;
d3 = b4*c2 + c4*b2;
d4 = c4 * b4 * (d2 * a4 + d1 * a2 + x2 * d3);
d5 = (a4 * d1 + b4 * c4 * x2);
z = a2*( d4.*x.*cos(x) +  a2*d1.*d5.*sin(x)) ./ d5.^2;
end

function z = u(p) 
% exact solution
p = project(p);
z =sin(p(:,1));
end

function z = gra_u(p) 

p = project(p);
x = p(:,1); y = p(:,2); z = p(:,3);
x2 = x.^2; y2 = y.^2; z2 = z.^2;

d1 = x2/a4 + y2/b4 + z2/c4 ;

z1 = (1 - x2 ./(a4*d1));
z2 = -x.*y./(a2*b2*d1);
z3 = -x.*z./(a2*c2*d1);
z = cos(x)*ones(1,3).*[z1,z2,z3];
end

end




