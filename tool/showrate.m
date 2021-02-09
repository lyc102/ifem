function r = showrate(N,err,varargin)
%% SHOWRATE rate of an err sequence
%
%  r = showrate(N,err) finds the number r such that err = N^r and plots the
%  err vs N in loglog scale.
% 
%  r = showrate(N,err,k) finds the number r such that err(k:end)=N(k:end)^r.
%  
%  r = showrate(N,err,k,str) the legend uses str as the name of err
%
%  r = showrate(N,err,LineSpec,values) 
%  The function accepts standard LineSpec line properties. 
%  For more example, see https://www.mathworks.com/help/matlab/ref/linespec.html
%  
%  r = showrate(N,err,'r') Plot the error curve in red. 
%  r = showrate(N,err,'r-d') Plot the error curve in red dashed lines. 
%  r = showrate(N,err,'linewidth',2) Plot the error curve in thick lines. 
% 
%  See also showrate2, showresult, showmesh, showsolution
%
%  Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~isempty(varargin)
    if isnumeric(varargin{1}); k = varargin{1}; end
    for i=1:length(varargin)
        if ischar(varargin{i})
            if any(regexp(varargin{i},'[-.*^]'))
                opt{1} = varargin{i};
            end
        end
    end
    
    if length(varargin)>=2 && ~isnumeric(varargin{1})
        if ~any(regexp(varargin{1},'[-.*^]'))
            opt = varargin(1:end);
        end
    elseif length(varargin)>=3 && isnumeric(varargin{1})
        if ~any(regexp(varargin{2},'[-.*^]'))
            opt = varargin(2:end);
        end
    elseif length(varargin)<2 && ischar(varargin{1})
        if ~any(regexp(varargin{i},'[-.*^]'))
            str = varargin{1};
        end
    elseif length(varargin) == 2 && isnumeric(varargin{1})
        if ~any(regexp(varargin{i},'[-.*^]'))
            str = varargin{2};
        end
    end
end
if ~exist('k', 'var'); k = 1; end
err(abs(err) <= 1e-15) = 1e-15; % Prevent the case err = 0, log(err) = -Inf.
p = polyfit(log(N(k:end)),log(err(k:end)),1);
r = single(p(1));
s = 0.8*err(k)/N(k)^r;
if exist('opt','var')
    h = loglog(N,err,opt{:});
    set(h,'linewidth',2);
    hold on
    h = loglog(N,s*N.^r);
    set(h,'linewidth',1,'linestyle','-.','marker','none');    
else
    loglog(N,err,'-*','linewidth',2);
    hold on
    loglog(N,s*N.^r,'-.','linewidth',1)
end
axis tight;
xlabel('Number of unknowns'); ylabel('Error');
title(['Rate of convergence is CN^{' num2str(r,2) '}'],'FontSize', 14);
if exist('str','var')
    h_legend = legend(str,['CN^{' num2str(r,2) '}'],'LOCATION','Best');
    set(h_legend,'FontSize', 14);
else 
    h_legend = legend('Error',['CN^{' num2str(r,2) '}'],'LOCATION','Best');
    set(h_legend,'FontSize', 14);
end