function [gx, gy, gz] = gaussPtetra(X1, X2, X3, X4, ng)

%% USAGE: generate Gaussian nodes on tetra mesh
%
% INPUTS:
% p --- nt-by-3 vector
% t --- nt-by-4 vector
% ng --- number of Gaussian nodes in an triangle
%            possible values = 1, 4, 5, 11, 15.
%
% OUTPUTS:
% gx --- nt-by-ng vector: stores x coordinates of Gaussian nodes
% gy --- nt-by-ng vector: stores y coordinates of Gaussian nodes
% gz --- nt-by-ng vector: stores z coordinates of Gaussian nodes
%
% Last Modified: 07/02/2020 by Xu Zhang
%%

G = zeros(size(X1,1),3*ng); i = 1;

if ng == 1 % accurate up to polynomial degree p = 1
    
    w1 = 1/4; % point at centroid
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X2+X3+X4);
    
elseif ng == 4 % accurate up to polynomial degree p = 2
    
    w1 = 0.585410196624969; w2 = 0.138196601125011; % median line    
    G(:,[i,i+ng,i+2*ng]) = w1*X1 + w2*(X2+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X2 + w2*(X1+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X3 + w2*(X1+X2+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X4 + w2*(X1+X2+X3);  
    
elseif ng == 5 % accurate up to polynomial degree p = 3
    
    w1 = 1/4; % point at centroid
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X2+X3+X4); i = i+1;
    
    w1 = 1/2; w2 = 1/6;% median line
    G(:,[i,i+ng,i+2*ng]) = w1*X1 + w2*(X2+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X2 + w2*(X1+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X3 + w2*(X1+X2+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X4 + w2*(X1+X2+X3);
    
elseif ng == 11 % accurate up to polynomial degree p = 4
    
    w1 = 1/4; % point at centroid
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X2+X3+X4); i = i+1;
    
    w1 = 0.785714285714286; w2 = 0.071428571428571; % median line
    G(:,[i,i+ng,i+2*ng]) = w1*X1 + w2*(X2+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X2 + w2*(X1+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X3 + w2*(X1+X2+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X4 + w2*(X1+X2+X3); i = i+1;
    
    w1 = 0.399403576166799; w2 = 0.100596423833201; % line opposite sides
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X2) + w2*(X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X3) + w2*(X2+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X4) + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X2+X3) + w2*(X1+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X2+X4) + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X3+X4) + w2*(X1+X2); 
    
elseif ng == 15 % accurate up to polynomial degree p = 5
    
    w1 = 1/4; % point at centroid
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X2+X3+X4); i = i+1;
    
    w1 = 0; w2 = 1/3; % median line
    G(:,[i,i+ng,i+2*ng]) = w1*X1 + w2*(X2+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X2 + w2*(X1+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X3 + w2*(X1+X2+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X4 + w2*(X1+X2+X3); i = i+1; 
    
    w1 = 0.727272727272727; w2 = 0.090909090909091; % median line
    G(:,[i,i+ng,i+2*ng]) = w1*X1 + w2*(X2+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X2 + w2*(X1+X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X3 + w2*(X1+X2+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*X4 + w2*(X1+X2+X3); i = i+1; 
    
    w1 = 0.066550153573664; w2 = 0.433449846426336; % line opposite sides
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X2) + w2*(X3+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X3) + w2*(X2+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X1+X4) + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X2+X3) + w2*(X1+X4); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X2+X4) + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng,i+2*ng]) = w1*(X3+X4) + w2*(X1+X2); 
    
end
gx = G(:,1:ng); gy = G(:,ng+1:2*ng); gz = G(:,2*ng+1:3*ng); 