function output = array2cell(x, rows)
%% array2cell: produces cell array using the cut 
%  A modification of mat2cell to output row-wise in each element
%  output cell array itself is column
if ~isequal(sum(rows),size(x,1))
   error('array2cell:VectorSumMismatch',...
       'split size is %s, but array size is %s.',num2str(sum(rows)), num2str(size(x,1)))
end
output = mat2cell(x',1,rows')';
end
