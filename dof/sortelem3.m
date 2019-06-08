function [elem,bdFlag] = sortelem3(elem,bdFlag)
%% SORTELEM3 sort elem in ascend ordering
%
% [elem,bdFlag] = sortelem3(elem,bdFlag) sorts the elem such that
% elem(t,1)< elem(t,2)< elem(t,3)<elem(t,4). A simple sort(elem,2) cannot
% sort bdFlag.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Step 1: elem(:,4) is the biggest one
[tempvar,idx] = max(elem,[],2);  %#ok<*ASGLU>
elem(idx==1,1:4) = elem(idx==1,[2 4 3 1]);
elem(idx==2,1:4) = elem(idx==2,[3 4 1 2]);
elem(idx==3,1:4) = elem(idx==3,[4 2 1 3]);
if exist('bdFlag','var')
    bdFlag(idx==1,1:4) = bdFlag(idx==1,[2 4 3 1]);
    bdFlag(idx==2,1:4) = bdFlag(idx==2,[3 4 1 2]);
    bdFlag(idx==3,1:4) = bdFlag(idx==3,[4 2 1 3]);
end
%% Step 2: elem(:,1) is the smallest one
[tempvar,idx] = min(elem(:,1:3),[],2);
% elem(idx==1,1:3) = elem(idx==1,[1 2 3]);
elem(idx==2,1:3) = elem(idx==2,[2 3 1]);
elem(idx==3,1:3) = elem(idx==3,[3 1 2]);
if exist('bdFlag','var')
    bdFlag(idx==2,1:3) = bdFlag(idx==2,[2 3 1]);
    bdFlag(idx==3,1:3) = bdFlag(idx==3,[3 1 2]);
end

%% Step 3: sort elem(:,2)<elem(:3)
idx = (elem(:,3) < elem(:,2));
elem(idx,[2 3]) = elem(idx,[3 2]);
if exist('bdFlag','var')
    bdFlag(idx,[2 3]) = bdFlag(idx,[3 2]); 
end

%% Output
if ~exist('bdFlag','var')
    bdFlag = [];
end