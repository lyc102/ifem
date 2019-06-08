function [node,elem] = squaregradmesh(square,h,s)
%
%   [node,elem] = squaregradmesh([0,1,0,2],0.1,0.2);
%   showmesh(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

% Generate a uniform mesh first
[node,elem] = squaremesh(square,h);

%% Change coordinate of y
x0 = square(1); x1 = square(2); 
y0 = square(3); y1 = square(4);
Mx = (x1-x0)/h;
Tx = x0:h:x1; % uniform in x
if s == 0.5
    r = 1;
else
    r = 3/(2*s)+0.1;
end
My = (y1-y0)/h;
k = 0:My;
Ty = (k/My).^r*(y1-y0) + y0;
x = kron(Tx,ones(My+1,1));
y = kron(ones(1,Mx+1),Ty');
node = [x(:),y(:)];