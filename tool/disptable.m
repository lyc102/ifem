function disptable(colname,varargin)
%% DISPTABLE display data in table
%
% DISPTABLE({'Data 1', 'Data 2'},data1,format1,data2,format2) display data1
% and data2 in table using format1 and format2, respectively.  
%
% Example
%
%	load dispdata
%   colname = {'#nodes'     '|u_I-u_h|_1'      '||p-p_h||'};
%   disptable(colname,err.N,[],err.uH1,'%0.5e',err.pL2,'%0.5e');
%
% See also displaytable
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

displaytable(colname,varargin{:});