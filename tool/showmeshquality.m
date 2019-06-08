function [q,u,N,x] = showmeshquality(node,elem,fh,varargin)
%% SHOWMESHQUALITY show the qualify of the mesh
%
% showmeshquality(node,elem) display the qualify of each triangle using
% histogram count and print on the screen.
%
% showmeshquality(node,elem,fh) also compute the uniformity of the mesh
% size
%
% Example
%   load lakemesh
%   showmeshquality(node,elem);
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

q = simpqual(node,elem);
if nargin >2
    [c,r] = circumcenter(node,elem);
    hc = feval(fh,c,varargin{:});
    sz = r./hc;
    u = std(sz)/mean(sz);
else
    u = 0;
end
x = 0:0.02:1;
hist(q,x);
[N,x] = histc(q,x);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.55 0.55 0.55],'EdgeColor','k')
axis tight
minstr = num2str(min(q),'%.4f\n');
meanstr = num2str(mean(q),'%.4f\n');
uniformstr = num2str(100*u,'%2.2f%%\n');
k = max([length(minstr),length(meanstr),length(uniformstr)]);
minstr(k+1) = ' ';
meanstr(k+1) = ' ';
uniformstr(k+1) = ' ';
if nargin >2
    text(0.25,max(N)/2,['- Min        ' minstr;  ...
                        '- Mean       ' meanstr; ...
                        '- Uniformity ' uniformstr])
else
    text(0.25,max(N)/2,['- Min        ' minstr;  ...
                        '- Mean       ' meanstr])
end
fprintf(' - Min quality %.4f',min(q));
fprintf(' - Mean quality %.4f',mean(q));
if nargin > 2
    fprintf(' - Uniformity %.2f%%',100*u);
end
disp(' ')                