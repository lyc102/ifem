function ifem(topic)
%% IFEM Display HTML documentation in the MATLAB web browser.
%  
%     IFEM, by itself, displays contents of documentatoin of iFEM.
%  
%     IFEM FUNCTION displays the HTML documentation for the IFEM
%     function FUNCTION. 
%  
%     Examples:
%        ifem tutorial
%        ifem introduction
%        ifem Lshape
%        ifem crack
%        ifem mgdoc
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if nargin == 0
    currentFolder = pwd;
    ifempath = which('ifem');
    ifemdocpath = ifempath(1:35);
    cd(ifemdocpath)
    help Contents
    cd(currentFolder)
else
    topic = [topic '.html'];
    web(['file:///' which(topic)])
end
