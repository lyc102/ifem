function [node,elem] = regpolygon(n,h)

% Initial mesh
dtheta = 2*pi/n;
theta = (dtheta-pi)/2:dtheta:(3*pi-dtheta)/2;
node = [cos(theta'), sin(theta')];
node(end+1,:) = 0;
elem = delaunayn(node);
elem = fixorder(node,elem);

% Refinement
k = ceil(1/h);
for i=1:k
    [node,elem] = uniformrefine(node,elem);
end
