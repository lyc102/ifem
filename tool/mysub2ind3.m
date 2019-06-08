function idx = mysub2ind3(siz,i,j,k)
    nr = siz(1); nc = siz(2); nv = siz(3);
    if (max(j)>nc) || (max(i)>nr) || (max(k)>nv)
        error(message('MATLAB:mysub2ind:IndexOutOfRange'));
    end
    idx = (k-1)*nr*nc + (j-1)*nr+i;
end