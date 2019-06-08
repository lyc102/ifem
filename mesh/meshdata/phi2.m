function z = phi3(p) % level set function
x1 = acos(-1/4)/4;
x0 = pi/2 - x1 - sin(4*x1);
r0 = 0.60125;
r1 = 0.24012;
x = p(:,1); y = p(:,2);
 

theta = atan2(y,x);
isNeg = theta < 0;
theta(isNeg) = theta(isNeg) + 2*pi;

z = zeros(size(p,1),1);
rp = sqrt(x.^2 + y.^2);

isSingle0 = theta >= 0 & theta < x0; % [0,x0]
isSingle1 = theta > pi/2 - x0 & theta < pi/2 + x0;%[pi/2 - x0, pi/2 + x0]
isSingle2 = theta > pi - x0 & theta < pi + x0;% [pi - x0, pi + x0];
isSingle3 = theta > 3*pi/2 - x0 & theta < 3*pi/2 +x0;% [3*pi/2 - x0, 3*pi/2 + x0]
isSingle4 = theta > 2*pi - x0 & theta <= 2*pi;% [2*pi - x0, 2*pi]

isThree0 = theta >= x0 & theta <= pi/2 - x0;% [x0,x1], [x1,pi/2-x1],[pi/2 - x1,pi/2-x0]
isThree1 = theta >= pi/2 + x0 & theta <= pi - x0; % [pi/2 + x0, pi/2 + x1], [ pi/2 + x1, pi - x1],[pi - x1, pi - x0]
isThree2 = theta >= pi + x0 & theta <= 3*pi/2 - x0;% [pi + x0, pi + x1],[pi + x1, 3*pi/2 - x1],[3*pi/2 - x1,3*pi/2 - x0]
isThree3 = theta >= 3*pi/2 + x0 & theta <= 2*pi - x0;% [3*pi/2 + x0, 3*pi/2 + x1], [ 3*pi/2 + x1, 2*pi - x1],[2*pi - x1, 2*pi - x0]

isSingle = isSingle0 | isSingle1 | isSingle2 | isSingle3 | isSingle4;
theta1 = theta(isSingle);
t0 = zeros(size(p,1),1);
t0(isSingle0) = x0/2;
t0(isSingle1) = pi/2;
t0(isSingle2) = pi;
t0(isSingle3) = 3*pi/2;
t0(isSingle4) = 2*pi - x0/2;
if ~isempty(theta1)
   t = getT(theta1,t0(isSingle));
   r = r0 + r1*cos(4*t + pi/2);
   z(isSingle) = rp(isSingle).^2 - r.^2;
end

isThree = isThree0 | isThree1 | isThree2 | isThree3;
theta1 = theta(isThree);
z1 = zeros(size(theta1,1),1);
theta1 = [theta1, theta1,theta1];
t0(isThree0,1) = (x0 + x1)/2; t0(isThree0,2) = pi/4; t0(isThree0,3) = pi/2 - (x0 + x1)/2;
t0(isThree1,1) = pi/2 + (x0 + x1)/2; t0(isThree1,2) = 3*pi/4; t0(isThree1,3) = pi - (x0 + x1)/2;
t0(isThree2,1) = pi + (x0 + x1)/2; t0(isThree2,2) = 5*pi/4; t0(isThree2,3) = 3*pi/2 - (x0 + x1)/2;
t0(isThree3,1) = 3*pi/2 + (x0 + x1)/2; t0(isThree3,2) = 7*pi/4; t0(isThree3,3) = 2*pi - (x0 + x1)/2;
if any(isThree)
   t = getT(theta1,t0(isThree,:));
   r = r0 + r1*cos(4*t + pi/2); 
   rt = rp(isThree);
   flag1 = rt < (r(:,1) + r(:,2))/2;
   flag2 = rt >= (r(:,1) + r(:,2))/2 & rt < (r(:,2) + r(:,3))/2;
   flag3 = rt >= (r(:,2) + r(:,3))/2 ;
   if any(flag1)
       z1(flag1) = rt(flag1).^2 - r(flag1,1).^2;
   end
   if any(flag2)
       z1(flag2) = r(flag2,2).^2 - rt(flag2).^2;
   end
   if any(flag3)
       z1(flag3) = rt(flag3).^2- r(flag3,3).^2;
   end
   z(isThree) = z1;
end

tt = [x1, 2*pi - x1, pi/2 - x1, pi/2 + x1, pi - x1, pi + x1, 3*pi/2 - x1, 3*pi/2 + x1];
rt = r0 + r1 * cos(4*tt + pi/2);
xt = rt.*cos( tt + sin(4*tt));
yt = rt.*sin( tt + sin(4*tt));
rt = zeros(length(x),8);
for i = 1:8
    rt(:,i) = sqrt((x - xt(i)).^2 + (y - yt(i)).^2);    
end
rt = min(rt,[],2);

z = msign(z).*min([abs(z),rt],[],2);
end


function t = getT(theta,t0)

f = t0 + sin(4*t0) - theta;
fprime = 1 + 4*cos(4*t0);
t = t0 - f./fprime;
err = sqrt(sum(sum(f.^2)));
while err > sqrt(eps)
    t0 = t;
    f = t0 + sin(4*t0) - theta;
    fprime = 1+4*cos(4*t0);
    t = t0 - f./fprime;
    err = sqrt(sum(sum(f.^2)));
end
end
