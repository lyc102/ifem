function [elem,bdFlag] = sortelem(elem,bdFlag)
%% SORTELEM sort elem in ascend ordering
%
% [elem,bdFlag] = SORTELEM(elem,bdFlag) sorts the elem such that
% elem(t,1)< elem(t,2)< elem(t,3). A simple sort(elem,2) cannot
% switch bdFlag.
%
% See also  sortelem3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Dimension check
if size(elem,2) >= 4 % 3D case
    [elem,bdFlag] = sortelem3(elem,bdFlag);
end

%% Step 1: make elem(:,3) to the largest one
[tempvar,idx] = max(elem,[],2); %#ok<ASGLU>
elem(idx==1,1:3) = elem(idx==1,[2 3 1]);
elem(idx==2,1:3) = elem(idx==2,[3 1 2]);
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlag(idx==1,1:3) = bdFlag(idx==1,[2 3 1]);
    bdFlag(idx==2,1:3) = bdFlag(idx==2,[3 1 2]);
end

%% Step 2: swtich the first two vertices such that elem(:,1)<elem(:2)
idx = (elem(:,2) < elem(:,1));
elem(idx,[1 2]) = elem(idx,[2 1]);
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlag(idx,[1 2]) = bdFlag(idx,[2 1]); 
else
    bdFlag = [];
end