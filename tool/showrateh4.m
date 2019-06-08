function showrateh4(h1,err1,k1,opt1,str1,h2,err2,k2,opt2,str2,...
                    h3,err3,k3,opt3,str3,h4,err4,k4,opt4,str4)
%% SHOWRATEH4 rate of two error sequences
%
% showrateh4(N1,err1,k1,opt1,str1,N2,err2,k2,opt2,str2,N3,err3,k3,opt3,str3,err4,k4,opt4,str4)
% plots the err1 vs N1, err2 vs N2, err3 vs N3, and err4 vs N4 in the
% loglog scale. Additional input
%
%   - k1, k2, k3, k4: specify the starting indices; see showrate
%   - opt1, opt2, opt3, opt4: the line color and style 
%   - str1, str2, str3, str4: strings used in legend
%
% Example
%
%
% See also showrate, showresult, showmesh, showsolution
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N1 = 1./h1; N2 = 1./h2; N3 =1./h3; N4 =1./h4;
if (nargin<=2) 
    k1 = 1; opt1 = '-*';
end
r1 = showrate(N1,err1,k1,opt1);
hold on
r2 = showrate(N2,err2,k2,opt2);
r3 = showrate(N3,err3,k3,opt3);
r4 = showrate(N4,err4,k4,opt4);
title(['Rate of convergence is Ch^{' num2str(r4,2) '}'],'FontSize', 14);
h_legend = legend(str1,['C_1h^{' num2str(-r1) '}'],...
                  str2,['C_2h^{' num2str(-r2) '}'],...
                  str3,['C_3h^{' num2str(-r3) '}'],...
                  str4,['C_4h^{' num2str(-r4) '}'],'LOCATION','Best');
set(h_legend,'FontSize',12);
xlabel('log(1/h)');
