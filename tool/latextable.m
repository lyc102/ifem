function latextable(A,filename)
%% LATEXTABLE output the matrix in table formate in latex.
%
% latextable(A,'A.tex');
% latextable(A);
% 
% If the file name is not given, the default one is 'table.tex'.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

m = size(A,1); n = size(A,2);
if (nargin <=2), filename='table.tex'; end %default file name
fid = fopen(filename,'wt');
fprintf(fid,'\\begin{tabular}{c} \n');
for i=1:m
    if i<m
        fprintf(fid,'%d \\\\ \n',i);
    else
        fprintf(fid,'%d \n',i);
    end
end
fprintf(fid,'\\end{tabular} \n');
fprintf(fid,'\\begin{tabular}{');
for j=1:n
    fprintf(fid,'|c');
end
fprintf(fid,'|} \n \\hline \n');
for i=1:m
    for j=1:n
        fprintf(fid,'%g',A(i,j));
        if j<n
            fprintf(fid,' & ');
        end
    end
    fprintf(fid,'\\\\ \n \\hline \n');
end
fprintf(fid,'\\end{tabular} \n \\\\');
fprintf(fid,'\\smallskip \n');
fprintf(fid,'\\begin{tabular}{');
for j=1:n
    fprintf(fid,'c');
end
fprintf(fid,'} \n \\quad \\;');
for j=1:n
    fprintf(fid,'%d', j);
    if j<n
        fprintf(fid,' & ');
    end
end
fprintf(fid,'\n \\end{tabular} \n');
fclose(fid);
%% TODO: add latexerrtable