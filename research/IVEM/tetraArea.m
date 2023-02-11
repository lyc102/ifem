function A = tetraArea(X1,X2,X3,X4)

%% Usage: calculate the area of each tetrahedral element
% INPUTS:
% X1 = nt-by-3 vector, stores [x,y,z] coordinate of 1st node of each cell 
% X2 = nt-by-3 vector, stores [x,y,z] coordinate of 2nd node of each cell 
% X3 = nt-by-3 vector, stores [x,y,z] coordinate of 3rd node of each cell 
% X4 = nt-by-3 vector, stores [x,y,z] coordinate of 4th node of each cell 

% OUTPUTS:
% A --- area of the triangle formed above

% Last Modified: 08/07/2020 by Xu Zhang
%%
a = X1 - X4; b = X2 - X4; c = X3 - X4;
A = (1/6)*abs(dot(a',cross(b',c')))';
