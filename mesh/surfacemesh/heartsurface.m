function surfacedata = heartsurface

surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'Hessian', @Hessian, ...
    'tanproperator', @tanproperator, 'initmesh',@initmesh, 'meancurvature', @meancurvature);


function [node,elem] = initmesh
M = load('heart2697.mat');
node = M.node;
elem = M.elem;

function z = phi(p)
% level set function
z =  (p(:,1) - p(:,3) .^ 2) .^ 2 + p(:,2) .^ 2 + p(:,3) .^ 2 - 1;


function n = unitoutnormal(p)
% 
s = 1 + 4 .* p(:,3) .^ 6 + (-8 .* p(:,1) + 4) .* p(:,3) .^ 4 + (-4 .* p(:,1) + 4 .* p(:,1) .^ 2) .* p(:,3) .^ 2;
n = [(s .^ (-0.1e1 ./ 0.2e1) .* (2 .* p(:,1) - 2 .* p(:,3) .^ 2)) ./ 0.2e1, s .^ (-0.1e1 ./ 0.2e1) .* p(:,2),...
    (s .^ (-0.1e1 ./ 0.2e1) .* (-4 .* (p(:,1) - p(:,3) .^ 2) .* p(:,3) + 2 .* p(:,3))) ./ 0.2e1];

function n = gradient(p)
% 
n = [2 .* p(:,1) - 2 .* p(:,3) .^ 2, 2 .* p(:,2), -4 .* (p(:,1) - p(:,3) .^ 2) .* p(:,3) + 2 .* p(:,3)];

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

function H = Hessian(p)

H = zeros(3,3,size(p,1));
s = 1 + 4 .* p(:,3) .^ 6 + (-8 .* p(:,1) + 4) .* p(:,3) .^ 4 + (-4 .* p(:,1) + 4 .* p(:,1) .^ 2) .* p(:,3) .^ 2;
d11 = (sqrt(s) - (p(:,1) .^ 2) + (4 .* p(:,1) .* p(:,3) .^ 2) - (4 .* p(:,1) .^ 2 .* p(:,3) .^ 2) + (8 .* p(:,1) .* p(:,3) .^ 4) - (3 .* p(:,3) .^ 4) - (4 .* p(:,3) .^ 6)) .* s .^ (-0.3e1 ./ 0.2e1);

d12 = -(p(:,1) - p(:,3) .^ 2) .* s .^ (-0.3e1 ./ 0.2e1) .* p(:,2);

d13 = -p(:,3) .* (2 .* s - 6 .* p(:,1) .^ 2 + 16 .* p(:,1) .* p(:,3) .^ 2 + 4 .* p(:,1) .^ 3 - 20 .* p(:,1) .^ 2 .* p(:,3) .^ 2 + 28 .* p(:,1) .* p(:,3) .^ 4 + p(:,1) - 10 .* p(:,3) .^ 4 - 12 .* p(:,3) .^ 6 - p(:,3) .^ 2) .* s .^ (-0.3e1 ./ 0.2e1);

d21 = -(p(:,2) .* s .^ (-0.3e1 ./ 0.2e1) .* (2 .* p(:,1) - 6 .* p(:,3) .^ 2 + 8 .* p(:,1) .* p(:,3) .^ 2 - 8 .* p(:,3) .^ 4)) ./ 0.2e1;

d22 = (s - p(:,2) .^ 2) .* s .^ (-0.3e1 ./ 0.2e1);

d23 = -p(:,2) .* p(:,3) .* (-6 .* p(:,1) + 10 .* p(:,3) .^ 2 + 4 .* p(:,1) .^ 2 - 16 .* p(:,1) .* p(:,3) .^ 2 + 12 .* p(:,3) .^ 4 + 1) .* s .^ (-0.3e1 ./ 0.2e1);

d31 = -p(:,3) .* (2 .* s - 2 .* p(:,1) .^ 2 + 12 .* p(:,1) .* p(:,3) .^ 2 - 8 .* p(:,1) .^ 2 .* p(:,3) .^ 2 + 16 .* p(:,1) .* p(:,3) .^ 4 - 10 .* p(:,3) .^ 4 - 8 .* p(:,3) .^ 6 + p(:,1) - 3 .* p(:,3) .^ 2) .* s .^ (-0.3e1 ./ 0.2e1);

d32 = p(:,3) .* (2 .* p(:,1) - 2 .* p(:,3) .^ 2 - 1) .* s .^ (-0.3e1 ./ 0.2e1) .* p(:,2);

d33 = -(2 .* s .* p(:,1) - 6 .* p(:,3) .^ 2 .* s - s + 16 .* p(:,1) ...
        .^ 2 .* p(:,3) .^ 2 - 48 .* p(:,1) .* p(:,3) .^ 4 - 8 .* p(:,3) .^ 2 .* p(:,1) .^ 3 + 40 .* p(:,3) .^ 4 .* p(:,1) .^ 2 - 56 .* p(:,3) .^ 6 .* p(:,1) - 8 .* p(:,1) .* p(:,3) .^ 2 + 32 .* p(:,3) .^ 6 + 24 .* p(:,3) .^ 8 + 12 .* p(:,3) .^ 4 + p(:,3) .^ 2) .* s .^ (-0.3e1 ./ 0.2e1);

for i=1:size(p,1)
    H(:,:,i)=[d11(i),d21(i),d31(i);d12(i),d22(i),d32(i);d13(i),d23(i),d33(i)];
end

function mc = meancurvature(p)
x = p(:,1);
y = p(:,2);
z = p(:,3);
s1 = -3 * z.^ 6 + (7 * x - 3) .* z.^ 4 + (4 * x - 5 * x.^ 2 - 1 - 3 * y.^ 2) .* z.^ 2 + (x - 1) .* (x.^ 2 + y.^2);
s2 = (4 * z.^ 6 + (5 - 8 * x) .* z.^ 4 + (-6 * x + 4 * x.^ 2 + 1).* z.^ 2 + x.^ 2 + y.^ 2).^ (0.3e1 / 0.2e1);
mc = -s1./s2;

function z = tanproperator(p)

z = zeros(3,3,size(p,1));
n = unitoutnormal(p);
for i=1:size(p,1)
   z(:,:,i)=n(i,:)'*n(i,:); 
end


