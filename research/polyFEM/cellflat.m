function out = cellflat(celllist)
% CELLFLAT is a helper function to flatten nested cell arrays. 
% 
% CELLFLAT(celllist) searches every cell element in cellist and put them on
% the top most level. Therefore, CELLFLAT linearizes a cell array tree
% structure. 
%
% Example: cellflat({[1 2 3], [4 5 6],{[7 8 9 10],[11 12 13 14 15]},{'abc',{'defg','hijk'},'lmnop'}}) 
% 
% Output: 
%Columns 1 through 7
%     [1x3 double]    [1x3 double]    [1x4 double]    [1x5 double]    'abc'    'defg'    'hijk'
%   Column 8 
%     'lmnop'
%
% cellflat(({{1 {2 3}} 'z' {'y' 'x' 'w'} {4 @iscell 5} 6}) )
% Output: 
% [1]    [2]    [3]    'z'    'y'    'x'    'w'    [4]    @iscell    [5]    [6]
%
% Version: 1.0
% Author: Yung-Yeh Chang, Ph.D. (yungyeh@hotmail.com)
% Date: 12/31/2014
% Copyright 2015, Yung-Yeh Chang, Ph.D.
% See Also: cell
if ~iscell(celllist)
    error('CELLFLAT:ImproperInputAugument','Input argument must be a cell array');
end
out = {};
for idx_c = 1:numel(celllist)
    if iscell(celllist{idx_c})
        out = [out cellflat(celllist{idx_c})];
    else
        out = [out celllist(idx_c)];
    end
end