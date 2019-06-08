function includepath(pathstr)
%% INCLUDEPATH inlcudes pathes
% 
% includepath(path) includes all subdirectories under pathstr to search
% path
%
% Example
%
%   ifempath = '/Users/longchen/Math/Programming/iFEM';
%   addpath(ifempath); 
%   includepath(ifempath);
%
% See also: setpath, addpath
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if nargin<1
    addpath(genpath(pwd),'-begin');
else
    currentFolder = pwd;
    cd(pathstr);
    addpath(genpath(pwd),'-begin');
    cd(currentFolder);
end