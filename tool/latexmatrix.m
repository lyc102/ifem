function latexmatrix(A,filename)
%% LATEXMATRIX output the matrix in latex formate.
%
% latexmatrix(A,'matrixA.tex');
% latexmatrix(A);
% 
% If the file name is not given, the default one is 'matrix.tex'.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

m = size(A,1); n = size(A,2);
if (nargin <=2), filename='matrix.tex'; end %default file name
fid = fopen(filename,'wt');
fprintf(fid,'\\left [ \\; \n');
fprintf(fid,'\\begin{matrix} \n');
for i=1:m
    for j=1:n
        fprintf(fid,'%d',A(i,j));
        if j<n
            fprintf(fid,' & ');
        else
            fprintf(fid,' \\, ');
        end
    end
    if i<m
        fprintf(fid,'\\\\ \n');
    else
        fprintf(fid,'\n');
    end
end
fprintf(fid,'\\end{matrix} \n');
fprintf(fid,'\\right ]');
fclose(fid);