function idx = mysub2ind(siz,i,j)

nr = siz(1); nc = siz(2);
if (max(j)>nc) || (max(i)>nr)
    error(message('MATLAB:mysub2ind:IndexOutOfRange'));
end
idx = (j-1)*nr+i;
