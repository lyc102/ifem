function [sortA,i2,j] = myunique(A)
%% MYUNIQUE the same input and output as unique 
%
% Solve the change of unique command in different version of MATLAB
%
% See also: unique
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

matlabversion = version;
startIndex = regexp(version,'(')+2;
endIndex = regexp(version, ')')-2;

if str2double(matlabversion(startIndex:endIndex)) <= 2012
    [sortA, i2, j] = unique(A,'rows');
else
    [sortA, i2, j] = unique(A,'rows','legacy');
end
