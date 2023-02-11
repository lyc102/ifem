function [gx, gy] = gaussPtri(p, t, ng)

%% USAGE: generate Gaussian nodes on a triangle
%
% INPUTS:
% p --- nt-by-2 vector
% t --- nt-by-3 vector
% ng --- number of Gaussian nodes in an triangle
%            possible values = 1, 3, 4, 6, 7, 12, 13.
%
% OUTPUTS:
% gx --- nt-by-ng vector: stores x coordinates of Gaussian nodes
% gy --- nt-by-ng vector: stores y coordinates of Gaussian nodes

% Last Modified: 07/05/2020 by Xu Zhang
%%
X1 = p(t(:,1),:);  
X2 = p(t(:,2),:);  
X3 = p(t(:,3),:); 
G = zeros(size(t,1),2*ng); i = 1;

if ng == 1 % accurate up to polynomial degree p = 1
    
    w1 = 1/3; % point at centroid
    G(:,[i,i+ng]) = w1*(X1+X2+X3);
    
elseif ng == 3 % accurate up to polynomial degree p = 2
    
    w1 = 2/3; w2 = 1/6; % median line
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); 
    
elseif ng == 4 % accurate up to polynomial degree p = 3
    
    w1 = 1/3; % point at centroid
    G(:,[i,i+ng]) = w1*(X1+X2+X3); i = i+1;
    
    w1 = 0.6; w2 = 0.2; % median line
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2);  
    
elseif ng == 6 % accurate up to polynomial degree p = 4
    
    w1 = 0.816847572980458; w2 = 0.091576213509771; 
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); i = i+1;
    
    w1 = 0.108103018168070; w2 = 0.445948490915965; 
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); 
    
elseif ng == 7 % accurate up to polynomial degree p = 5
    
    w1 = 1/3; % point at centroid
    G(:,[i,i+ng]) = w1*(X1+X2+X3); i = i+1;
    
    w1 = 0.797426985353087; w2 = 0.101286507323456;
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); i = i+1;
    
    w1 = 0.059715871789770; w2 = 0.470142064105115; 
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); 
    
elseif ng == 12 % accurate up to polynomial degree p = 6
    
    w1 = 0.873821971016996; w2 = 0.063089014491502; 
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); i = i+1;
    
    w1 = 0.501426509658180; w2 = 0.249286745170910; 
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); i = i+1;
   
    w1 = 0.636502499121399; w2 = 0.310352451033784; w3 = 0.053145049844817;
    G(:,[i,i+ng]) = w1*X1 + w2*X2 + w3*X3; i = i+1;
    G(:,[i,i+ng]) = w1*X1 + w3*X2 + w2*X3; i = i+1;
    G(:,[i,i+ng]) = w2*X1 + w1*X2 + w3*X3; i = i+1;
    G(:,[i,i+ng]) = w2*X1 + w3*X2 + w1*X3; i = i+1;
    G(:,[i,i+ng]) = w3*X1 + w1*X2 + w2*X3; i = i+1;
    G(:,[i,i+ng]) = w3*X1 + w2*X2 + w1*X3; 

elseif ng == 13 % accurate up to polynomial degree p = 7
    
    w1 = 1/3; % point at centroid
    G(:,[i,i+ng]) = w1*(X1+X2+X3); i = i+1;
    
    w1 = 0.479308067841920; w2 = 0.260345966079040; 
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); i = i+1;
   
    w1 = 0.869739794195568; w2 = 0.065130102902216; 
    G(:,[i,i+ng]) = w1*X1 + w2*(X2+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X2 + w2*(X1+X3); i = i+1;
    G(:,[i,i+ng]) = w1*X3 + w2*(X1+X2); i = i+1;
   
    w1 = 0.638444188569810; w2 = 0.312865496004874; w3 = 0.048690315425316;
    G(:,[i,i+ng]) = w1*X1 + w2*X2 + w3*X3; i = i+1;
    G(:,[i,i+ng]) = w1*X1 + w3*X2 + w2*X3; i = i+1;
    G(:,[i,i+ng]) = w2*X1 + w1*X2 + w3*X3; i = i+1;
    G(:,[i,i+ng]) = w2*X1 + w3*X2 + w1*X3; i = i+1;
    G(:,[i,i+ng]) = w3*X1 + w1*X2 + w2*X3; i = i+1;
    G(:,[i,i+ng]) = w3*X1 + w2*X2 + w1*X3; 
    
end
gx = G(:,1:ng); gy = G(:,ng+1:2*ng); 