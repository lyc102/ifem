function showrateh2(h1,err1,k1,opt1,str1,h2,err2,k2,opt2,str2)
%% SHOWRATEH2 rate of two error sequences
%
% showrateh2(N1,err1,k1,opt1,str1,N2,err2,k2,opt2,str2)
% plots the err1 vs N1, err2 vs N2, and err3 vs N3 in the loglog scale.
% Additional input
%
%   - k1, k2: specify the starting indices; see showrate
%   - opt1, opt2: the line color and style 
%   - str1, str2: strings used in legend
%
% Example
%
% showrateh2(h,energyErr,1,'r-+','||u-u_h||_A',...
%            h,L2Err,1,'b-+','||u-u_h||');
%
% See also showrateh, showresult, showmesh, showsolution
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N1 = 1./h1; N2 = 1./h2;
if (nargin<=2) 
    k1 = 1; opt1 = '-*';
end
r1 = -showrate(N1,err1,k1,opt1);
hold on
r2 = -showrate(N2,err2,k2,opt2);
title(['Rate of convergence is h^{' num2str(r2,2) '}'],'FontSize', 14);
h_legend = legend(str1,['C_1 h^{' num2str(r1,2) '}'],...
                  str2,['C_2 h^{' num2str(r2,2) '}'],'LOCATION','Best');
set(h_legend,'FontSize',12);
xlabel('log(1/h)');
