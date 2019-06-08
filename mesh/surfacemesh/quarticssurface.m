function surfacedata = quarticssurface

surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'initmesh',@initmesh,...
    'meancurvature', @meancurvature, 'u', @u, 'f', @f, 'gra_u',@gra_u);


function [node, elem] = initmesh
M = load('meshdata/quartics4528.mat');
node = M.node;
elem = M.elem;
end


function z = phi(p)
% level set function
x = p(:,1); y = p(:,2); z = p(:,3);
z = (x.^2 - 1).^2 + (y.^2 - 1).^2 + (z.^2 - 1).^2 - 1.05;
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
n = [4*(x.^2 - 1).*x,4*(y.^2 - 1).*y, 4*(z.^2 - 1).*z];
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
s = (x.^2-1).^2.*x.^2 + (y.^2-1).^2.*y.^2 + (z.^2-1).^2.*z.^2;
s1 =(3 .* z .^ 2 - 2 + 3 .* y .^ 2) .* x .^ 6 + (-6 .* y .^ 2 - 6 .* z .^ 2 + 4) .* x .^ 4 ...
    + (-2 + 6 .* y .^ 2 + 3 .* z .^ 6 - 6 .* z .^ 4 - 6 .* y .^ 4 + 6 .* z .^ 2 + 3 .* y .^ 6) .* x .^ 2 ...
+ (3 .* z .^ 2 - 2) .* y .^ 6 + (4 - 6 .* z .^ 2) .* y .^ 4 + ...
    (-2 + 6 .* z .^ 2 + 3 .* z .^ 6 - 6 .* z .^ 4) .* y .^ 2 - 2 .* z .^ 6 - 2 .* z .^ 2 + 4 .* z .^ 4;

mc = 0.5 * s1 ./ (s.^(3/2));
end


function z = f(p)
p = project(p);
x = p(:,1); y = p(:,2); z = p(:,3);
x2 = x.^2; y2 = y.^2; z2 = z.^2;
x4 = x2.^2; y4 = y2.^2; z4 = z2.^2;

x31 = x2.*x - x;
y31 = y2.*y - y;
z31 = z2.*z - z;

d1 = x31.^2 + y31.^2 + z31.^2;
d2 = x31 + y31 + z31;
t = exp(x+y+z);

d3 = (z2 + y2 - 2/3).*x4.*(x2 - 2);
d4 = (z4.*(z2-2) + 2*z2 + y4.*(y2 -2) + 2*y2 -2/3).*x2;
d5 = (z2 - 2/3).*y4.*(y2 - 2);
d6 = (z4.*(z2 -2) + 2*z2 - 2/3).*y2;

z = -2*t.*((-(y31+z31).*x31 - z31.*y31 + d1).*d1 -1.5*(d3 + d4 + d5 + d6 - 2/3*z31.^2).*d2)./d1.^2;
end

function z = u(p) 
% exact solution
p = project(p);
z = exp(sum(p,2));
end

function z = gra_u(p) 

p = project(p);
x = p(:,1); y = p(:,2); z = p(:,3);
x2 = x.^2; y2 = y.^2; z2 = z.^2;


x31 = x2.*x - x;
y31 = y2.*y - y;
z31 = z2.*z - z;

d1 = x31.^2 + y31.^2 + z31.^2;
t = exp(x+y+z);
z1 = (y31+z31).*x31 - y31.^2 - z31.^2;
z2 = (x31+z31).*y31 - x31.^2 - z31.^2;
z3 = (x31+y31).*z31 - x31.^2 - y31.^2;

z = -(t./d1)*ones(1,3).*[z1,z2,z3];
end

end