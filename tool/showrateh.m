function r = showrateh(h,err,k,opt,str)
%% SHOWRATEH rate of an err sequence
%
%  r = SHOWRATEH(h,err) finds the number r such that err = N^r and plots the
%  err vs N in loglog scale.
% 
%  r = SHOWRATEH(h,err,k) finds the number r such that err(k:end)=N(k:end)^r.
%
%  The function accepts standard plotting properting. For example, r =
%  showrate(N,err,[],'r') will plot the error curve in red. 
%
% See also showrate2, showresult, showmesh, showsolution
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = 1./h;
if (nargin<=2) 
    k = 1; opt = '-*';
end
if (nargin<=3) 
    opt = '-*';
end
r = -showrate(N,err,k,opt);
if exist('str','var')
    h_legend = legend(str,['C_1h^{' num2str(r,2) '}'],'LOCATION','Best');
    set(h_legend,'FontSize', 14);
end
xlabel('log(1/h)');
title(['Rate of convergence is Ch^{' num2str(r,2) '}'],'FontSize', 14);