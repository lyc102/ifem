function [gw, gx, gy, gz] = gaussPedge(X1, X2, ng)

if ng == 1
 
    gw_ref = 2; g(1) = 0;
 
elseif ng == 2
 
    gw_ref = [1,1];
    g(1) = -0.57735026918962576451;
    g(2) = 0.57735026918962576451;
 
elseif ng == 3
 
    gw_ref=[5/9 8/9 5/9];
    g(1) = -0.77459666924148337704;
    g(2) = 0;
    g(3) = 0.77459666924148337704;
 
 
elseif ng == 4
    
    gw_ref = [0.347854845137454, 0.652145154862546, ...
        0.652145154862546, 0.347854845137454];
    g(1) = -sqrt(3+2*sqrt(6/5))/sqrt(7);
    g(2) = -sqrt(3-2*sqrt(6/5))/sqrt(7);
    g(3) = -g(2);
    g(4) = -g(1);
 
elseif ng == 5
 
    gw_ref = zeros(1,n_gnodes);
    gw_ref(1) = 0.2369268850561890875;
    gw_ref(2) = 0.478628670499366468;
    gw_ref(3) = 0.56888888888888888889;
    gw_ref(4) = gw_ref(2);
    gw_ref(5) = gw_ref(1);
    g(1) = -0.9061798459386639928;
    g(2) = -0.53846931010568309104;
    g(3) = 0;
    g(4) = -g(2);
    g(5) = -g(1);
 
elseif n_gnodes == 6
    
    gw_ref = zeros(1,n_gnodes);
    gw_ref(1) = 0.17132449237917034504;
    gw_ref(2) = 0.36076157304813860757;
    gw_ref(3) = 0.46791393457269104739;
    gw_ref(4) = gw_ref(3);
    gw_ref(5) = gw_ref(2);
    gw_ref(6) = gw_ref(1);
    g(1) = -0.93246951420315202781;
    g(2) = -0.66120938646626451366;
    g(3) = -0.23861918608319690863;
    g(4) = 0.23861918608319690863;
    g(5) = -g(2);
    g(6) = -g(1);
    
else
    
    disp('number of Gaussian nodes should be between 1 and 6.')
    
end

gx = zeros(size(X1,1),ng); gy = gx; gz = gx;
gw = (sum((X2-X1).^2,2).^(1/2)/2)*gw_ref;
for j = 1:ng
gj = (X2 - X1)/2*g(j) + (X2 + X1)/2;
gx(:,j) = gj(:,1); gy(:,j) = gj(:,2); gz(:,j) = gj(:,3);
end
% if abs(p1(1) - p2(1)) < 1e-14
%     gx = p1(1)*ones(1,n_gnodes);
%     gy = 0.5*(p2(2)-p1(2))*g + 0.5*(p2(2)+p1(2));
% else
%     gx = 0.5*(p2(1)-p1(1))*g + 0.5*(p2(1)+p1(1));
%     gy = p1(2)+(p2(2)-p1(2))*(gx-p1(1))/(p2(1)-p1(1));
% end


