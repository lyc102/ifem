function [mina,maxa] = getangle(node,elem)
%% GETANGLE get the minimum and maximum angles of a triangulation 
%
%
% 
% Author: Huayi Wei <huayiwei1984@gmail.com>
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

v1 = node(elem(:,3),:) - node(elem(:,2),:);
v2 = node(elem(:,1),:) - node(elem(:,3),:);
v3 = node(elem(:,2),:) - node(elem(:,1),:);

s1 = sqrt(sum(v1.^2,2));
s2 = sqrt(sum(v2.^2,2));
s3 = sqrt(sum(v3.^2,2));

a1=acosd((s2.^2 + s3.^2 - s1.^2)./(2*s2.*s3));
a2=acosd((s3.^2 + s1.^2 - s2.^2)./(2*s3.*s1));
a3=acosd((s1.^2 + s2.^2 - s3.^2)./(2*s1.*s2));

mina = min(min([a1,a2,a3],[],2));
maxa = max(max([a1,a2,a3],[],2));


