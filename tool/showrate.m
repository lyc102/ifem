function r = showrate(N,err,k,opt,str)
%% SHOWRATE rate of an err sequence
%
%  r = showrate(N,err) finds the number r such that err = N^r and plots the
%  err vs N in loglog scale.
% 
%  r = showrate(N,err,k) finds the number r such that err(k:end)=N(k:end)^r.
%
%  The function accepts standard plotting properting. For example, r =
%  showrate(N,err,[],'r') will plot the error curve in red. 
%
% See also showrate2, showresult, showmesh, showsolution
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('k','var') || isempty(k)
    k = 1; 
end
err(err == 0) = 1e-16; % Prevent the case err = 0, log(err) = -Inf.
p = polyfit(log(N(k:end)),log(err(k:end)),1);
r = single(p(1));
s = 0.75*err(1)/N(1)^r;
if exist('opt','var')
    h = loglog(N,err,opt);
    set(h,'linewidth',2);
    hold on
    h = loglog(N,s*N.^r,opt);
    set(h,'linewidth',1,'linestyle','--','marker','none');    
else
    loglog(N,err,'-*','linewidth',2);
    hold on
    loglog(N,s*N.^r,'k-.','linewidth',1)
end
axis tight;
xlabel('Number of unknowns'); ylabel('Error');
title(['Rate of convergence is CN^{' num2str(r,2) '}'],'FontSize', 14);
if exist('str','var')
    h_legend = legend(str,['CN^{' num2str(r,2) '}'],'LOCATION','Best');
    set(h_legend,'FontSize', 14);
end