function hg(varargin)
%% HG  the interface of command hg in matlab
%
%   
% Example:
%     hg clone  
%     hg pull
%     hg ci -m "some comment"
%     hg merge
%
%   Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.
str = 'hg';
for i = 1:nargin
  str = [str,' ', varargin{i}];
end
system(str);