function q = quadedge(node,edge,func1,func2)
%q = \int_e f ds
NE = size(edge, 1);
b = [1/2; (1 + sqrt(3/5))/2; (1 - sqrt(3/5))/2];
w = [4/9; 5/18; 5/18];

ve = [node(edge(:,2),1) - node(edge(:,1),1) , node(edge(:,2),2) - node(edge(:,1),2)];
length_ve = sqrt(ve(:,1).^2 + ve(:,2).^2);

n1 = edge(:,1);	 n2 = edge(:,2);
x1 = node(n1,1); y1 = node(n1,2);   
x2 = node(n2,1); y2 = node(n2,2); 

xx = x1*b'+x2*(1 - b)';
yy = y1*b'+y2*(1 - b)';

if nargin==4
    func = @(p) func1(p).*func2(p);
else
    func = @(p) func1(p);
end

f = func([xx(:) yy(:)]);
f = reshape(f, [NE, 3]);

q = length_ve.*(w(1)*f(:,1)+w(2)*f(:,2)+w(3)*f(:,3));